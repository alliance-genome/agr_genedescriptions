import logging

from ontobio import Ontology

from genedescriptions.commons import DataType, Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from pipelines.alliance.ateam_api_helper import get_anatomy_ontologies_roots, get_ontology_node_children
from pipelines.alliance.ateam_db_helper import get_expression_annotations, get_ontology_pairs, get_gene_data

logger = logging.getLogger(__name__)


class AllianceDataManager(DataManager):
    """
    This class is a subclass of DataManager and is used to manage data for the Alliance.
    It provides methods to load and save data, as well as to perform various data processing tasks.
    """

    def __init__(self, config: GenedescConfigParser, alliance_release_version: str = None):
        self.config = config
        self.anatomy_ontologies_roots = None
        super().__init__()

    def load_annotations_from_persistent_store(self, associations_type: DataType, taxon_id: str, provider: str):
        if associations_type == DataType.GO:
            pass
        elif associations_type == DataType.EXPR:
            associations = []
            for gea in get_expression_annotations(taxon_id=taxon_id):
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

    def load_ontology_from_ateam_api(self, ontology_type: DataType, provider: str = None):
        curie_prefix = ""
        ontology = Ontology()
        if ontology_type == DataType.GO:
            pass
        elif ontology_type == DataType.EXPR:
            if provider == "WB":
                curie_prefix = "WBbt"
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
            visited_nodes = set()
            while nodes:
                logger.debug(f"Number of visited nodes: {str(len(visited_nodes))}")
                node = nodes.pop(0)
                if node["descendantCount"] > 0:
                    children = get_ontology_node_children(node_curie=node["curie"])
                    for child in children:
                        if child["curie"] not in visited_nodes:
                            self.add_node_to_ontobio_ontology_if_not_exists(
                                term_id=child["curie"],
                                term_label=child["name"],
                                term_type="anatomy",
                                is_obsolete=False,
                                ontology=ontology,
                                check_exists=False)
                            nodes.append(child)
                            visited_nodes.add(child["curie"])
                        ontology.add_parent(id=child["curie"], pid=node["curie"], relation="subClassOf")

    def load_ontology_from_persistent_store(self, ontology_type: DataType, provider: str = None):
        curie_prefix = ""
        ontology = Ontology()
        if ontology_type == DataType.GO:
            pass
        elif ontology_type == DataType.EXPR:
            if provider == "WB":
                curie_prefix = "WBbt"
            ontology_pairs = get_ontology_pairs(curie_prefix=curie_prefix)
            for onto_pair in ontology_pairs:
                self.add_node_to_ontobio_ontology_if_not_exists(
                    term_id=onto_pair["parent_curie"],
                    term_label=onto_pair["parent_name"],
                    term_type=onto_pair["parent_type"],
                    is_obsolete=onto_pair["parent_is_obsolete"],
                    ontology=ontology)
                self.add_node_to_ontobio_ontology_if_not_exists(
                    term_id=onto_pair["child_curie"],
                    term_label=onto_pair["child_name"],
                    term_type=onto_pair["child_type"],
                    is_obsolete=onto_pair["child_is_obsolete"],
                    ontology=ontology)
                ontology.add_parent(id=onto_pair["child_curie"], pid=onto_pair["parent_curie"],
                                    relation="subClassOf" if onto_pair["rel_type"] == "IS_A" else "BFO:0000050")
            self.set_ontology(ontology_type=ontology_type, ontology=ontology, config=self.config)

    def load_gene_data_from_persistent_store(self, provider: str):
        genes = get_gene_data(provider=provider)
        for gene in genes:
            self.gene_data[gene["gene_id"]] = Gene(gene["gene_id"], gene["gene_symbol"], False, False)
