import concurrent.futures
import logging

from ontobio import Ontology

from genedescriptions.commons import DataType, Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from pipelines.alliance.ateam_api_helper import get_anatomy_ontologies_roots, get_ontology_node_children, \
    get_expression_annotations_from_api, get_data_providers_from_api
from pipelines.alliance.ateam_db_helper import get_expression_annotations, get_ontology_pairs, get_gene_data, \
    get_data_providers

logger = logging.getLogger(__name__)


provider_to_expression_curie_prefix = {
    "WB": "WBbt",
    "ZFIN": "ZFA",
    "FB": "FBbt",
    "MGI": "EMAPA",
    "XBXL": "XAO",
    "XBXT": "XAO"
}


class AllianceDataManager(DataManager):
    """
    This class is a subclass of DataManager and is used to manage data for the Alliance.
    It provides methods to load and save data, as well as to perform various data processing tasks.
    """

    def __init__(self, config: GenedescConfigParser, alliance_release_version: str = None):
        self.config = config
        self.anatomy_ontologies_roots = None
        super().__init__()

    def load_annotations(self, associations_type: DataType, taxon_id: str, provider: str, source: str = "db"):
        if associations_type == DataType.GO:
            pass
        elif associations_type == DataType.EXPR:
            associations = []
            if source == "db":
                gea_list = get_expression_annotations(taxon_id=taxon_id)
            elif source == "api":
                gea_list = get_expression_annotations_from_api(data_provider=provider)
            else:
                raise ValueError("source must be either 'db' or 'api'")
            for gea in gea_list:
                associations.append(DataManager.create_annotation_record(
                    source_line="",
                    gene_id=gea["gene_id"],
                    gene_symbol=gea["gene_symbol"],
                    gene_type="gene",
                    taxon_id=taxon_id,
                    object_id=gea["anatomy_id"],
                    qualifiers=["Verified"],
                    aspect="A",
                    ecode="EXP",
                    references="",
                    prvdr=provider,
                    date=None))
            self.expression_associations = DataManager.create_annot_set_from_legacy_assocs(
                assocs=associations, ontology=self.expression_ontology)
            self.expression_associations = self.remove_blacklisted_annotations(
                association_set=self.expression_associations, ontology=self.expression_ontology,
                terms_blacklist=self.config.get_module_property(module=Module.EXPRESSION,
                                                                prop=ConfigModuleProperty.EXCLUDE_TERMS))

    @staticmethod
    def add_node_to_ontobio_ontology_if_not_exists(term_id, term_label, term_type, is_obsolete, ontology,
                                                   check_exists: bool = True):
        """Add Term to Ontobio Ontology If Not Exists."""
        if not check_exists or (not ontology.has_node(term_id) and term_label):
            if is_obsolete in ["true", "True"]:
                meta = {
                    "deprecated": True, "basicPropertyValues": [
                        {"pred": "OIO:hasOBONamespace", "val": term_type}]
                }
            else:
                meta = {
                    "basicPropertyValues": [
                        {"pred": "OIO:hasOBONamespace", "val": term_type}]
                }
            ontology.add_node(id=term_id, label=term_label, meta=meta)

    def load_ontology(self, ontology_type: DataType, provider: str = None, source: str = "db"):
        if source not in ["db", "api"]:
            raise ValueError("source must be either 'db' or 'api'")

        ontology = Ontology()
        curie_prefix = self._get_curie_prefix(ontology_type, provider)
        if not curie_prefix:
            return None

        if ontology_type == DataType.EXPR and source == "api":
            self._load_expression_ontology_from_api(ontology, curie_prefix, provider)
        elif source == "db":
            self._load_ontology_from_db(ontology, curie_prefix)

        self._add_artificial_nodes(ontology, ontology_type, provider)
        self.set_ontology(ontology_type=ontology_type, ontology=ontology, config=self.config)

    @staticmethod
    def _get_curie_prefix(ontology_type: DataType, provider: str) -> str:
        if ontology_type == DataType.GO:
            return "GO"
        elif ontology_type == DataType.EXPR:
            return provider_to_expression_curie_prefix.get(provider, "")
        return ""

    def _load_expression_ontology_from_api(self, ontology, curie_prefix: str, provider: str):
        if self.anatomy_ontologies_roots is None:
            self.anatomy_ontologies_roots = get_anatomy_ontologies_roots()
        roots = [root for root in self.anatomy_ontologies_roots if root["curie"].startswith(curie_prefix)]
        for root in roots:
            self.add_node_to_ontobio_ontology_if_not_exists(
                term_id=root["curie"],
                term_label=root["name"],
                term_type="anatomy",
                is_obsolete=False,
                ontology=ontology,
                check_exists=False)
        nodes = roots
        visited_nodes = set(root["curie"] for root in roots)

        def process_node(node):
            children = get_ontology_node_children(node_curie=node["curie"])
            new_children = []
            for child in children:
                if child["curie"] not in visited_nodes:
                    self.add_node_to_ontobio_ontology_if_not_exists(
                        term_id=child["curie"],
                        term_label=child["name"],
                        term_type="anatomy",
                        is_obsolete=False,
                        ontology=ontology,
                        check_exists=False)
                    new_children.append(child)
                    visited_nodes.add(child["curie"])
                ontology.add_parent(id=child["curie"], pid=node["curie"], relation="subClassOf")
            return new_children

        with concurrent.futures.ThreadPoolExecutor() as executor:
            while nodes:
                logger.debug(f"Number of visited nodes: {len(visited_nodes)}")
                future_to_node = {executor.submit(process_node, node): node for node in nodes}
                nodes = []
                for future in concurrent.futures.as_completed(future_to_node):
                    nodes.extend([node for node in future.result() if node["childCount"] > 0])

    def _load_ontology_from_db(self, ontology, curie_prefix: str):
        ontology_pairs = get_ontology_pairs(curie_prefix=curie_prefix)
        added_nodes = set()
        for onto_pair in ontology_pairs:
            if onto_pair["parent_curie"] not in added_nodes:
                self.add_node_to_ontobio_ontology_if_not_exists(
                    term_id=onto_pair["parent_curie"],
                    term_label=onto_pair["parent_name"],
                    term_type=onto_pair["parent_type"],
                    is_obsolete=onto_pair["parent_is_obsolete"],
                    ontology=ontology,
                    check_exists=False)
                added_nodes.add(onto_pair["parent_curie"])
            if onto_pair["child_curie"] not in added_nodes:
                self.add_node_to_ontobio_ontology_if_not_exists(
                    term_id=onto_pair["child_curie"],
                    term_label=onto_pair["child_name"],
                    term_type=onto_pair["child_type"],
                    is_obsolete=onto_pair["child_is_obsolete"],
                    ontology=ontology,
                    check_exists=False)
                added_nodes.add(onto_pair["child_curie"])
            ontology.add_parent(id=onto_pair["child_curie"], pid=onto_pair["parent_curie"],
                                relation="subClassOf" if onto_pair["rel_type"] == "IS_A" else "BFO:0000050")

    def _add_artificial_nodes(self, ontology, ontology_type: DataType, provider: str):
        if ontology_type == DataType.EXPR and provider == "MGI":
            self.add_node_to_ontobio_ontology_if_not_exists("EMAPA_ARTIFICIAL_NODE:99999",
                                                            "embryo",
                                                            "anatomical_structure",
                                                            False,
                                                            ontology)
            ontology.add_parent("EMAPA_ARTIFICIAL_NODE:99999", "EMAPA:0", relation="subClassOf")
            self.add_node_to_ontobio_ontology_if_not_exists("EMAPA_ARTIFICIAL_NODE:99998",
                                                            "head",
                                                            "anatomical_structure",
                                                            False,
                                                            ontology)
            ontology.add_parent("EMAPA_ARTIFICIAL_NODE:99998", "EMAPA:0", relation="subClassOf")
            self.add_node_to_ontobio_ontology_if_not_exists(
                "EMAPA_ARTIFICIAL_NODE:99997",
                "gland",
                "anatomical_structure",
                False,
                ontology)
            ontology.add_parent("EMAPA_ARTIFICIAL_NODE:99997", "EMAPA:0", relation="subClassOf")
        elif ontology_type == DataType.EXPR and provider == "FB":
            self.add_node_to_ontobio_ontology_if_not_exists(
                "FBbt_ARTIFICIAL_NODE:99999",
                "organism",
                "",
                False,
                ontology)
            ontology.add_parent("FBbt_ARTIFICIAL_NODE:99999",
                                "FBbt:10000000",
                                relation="subClassOf")

    def load_gene_data(self, species_taxon: str, source: str = "db"):
        if source not in ["db", "api"]:
            raise ValueError("source must be either 'db' or 'api'")
        if source == "db":
            genes = get_gene_data(species_taxon=species_taxon)
            self.gene_data = {}
            for gene in genes:
                self.gene_data[gene["gene_id"]] = Gene(gene["gene_id"], gene["gene_symbol"], False, False)
        elif source == "api":
            raise NotImplementedError("API loading is not implemented yet")

    @staticmethod
    def load_data_providers(source: str = "db"):
        if source not in ["db", "api"]:
            raise ValueError("source must be either 'db' or 'api'")
        if source == "api":
            logger.info("Loading data providers from API")
            return get_data_providers_from_api()
        elif source == "db":
            logger.info("Loading data providers from DB")
            return get_data_providers()
        return None
