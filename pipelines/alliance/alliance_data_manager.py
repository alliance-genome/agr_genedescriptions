import concurrent.futures
import logging
import os
import requests
import tempfile

from ontobio import Ontology

from genedescriptions.commons import DataType, Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from pipelines.alliance.ateam_api_helper import get_ontology_roots, get_ontology_node_children, \
    get_expression_annotations_from_api, get_data_providers_from_api
from pipelines.alliance.ateam_db_helper import get_expression_annotations, get_ontology_pairs, get_gene_data, \
    get_data_providers, get_disease_annotations, get_best_human_orthologs_for_taxon

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

    def _load_go_annotations(self, provider: str):
        if provider in ["XBXT", "XBXL"]:
            provider = "XB"
        release_version = os.environ.get("ALLIANCE_RELEASE_VERSION")
        if not release_version:
            raise RuntimeError("ALLIANCE_RELEASE_VERSION not set in environment")
        fms_url = f"https://fms.alliancegenome.org/api/snapshot/release/{release_version}"
        response = requests.get(fms_url)
        if response.status_code != 200:
            raise RuntimeError(f"Failed to fetch FMS API: {response.status_code}")
        snapshot = response.json().get("snapShot", {})
        data_files = snapshot.get("dataFiles", [])
        gaf_file = None
        for f in data_files:
            if f.get("dataType", {}).get("name") == "GAF" and f.get("dataSubType", {}).get("name") == provider:
                gaf_file = f
                break
        if not gaf_file:
            raise RuntimeError(f"No GAF file found for provider {provider} in release {release_version}")
        gaf_url = gaf_file.get("s3Url") or gaf_file.get("stableURL")
        if not gaf_url:
            raise RuntimeError(f"No download URL found for GAF file for provider {provider}")
        logger.info(f"GAF file for provider {provider}: {gaf_url}")
        # Download the GAF file to a temp location
        with tempfile.NamedTemporaryFile(delete=False, suffix=".gaf.gz") as tmp_gaf:
            gaf_response = requests.get(gaf_url)
            tmp_gaf.write(gaf_response.content)
            tmp_gaf_path = tmp_gaf.name
        # Use DataManager's load_associations_from_file to load the GAF file
        self.load_associations_from_file(
            associations_type=DataType.GO,
            associations_url=gaf_url,
            associations_cache_path=tmp_gaf_path,
            config=self.config
        )
        os.remove(tmp_gaf_path)
        renamed_associations = []
        for subj_associations in self.go_associations.associations_by_subj.values():
            for association in subj_associations:
                parts = association["subject"]["id"].split(":")
                if len(parts) > 2 and parts[0] == parts[1]:
                    association["subject"]["id"] = f"{parts[0]}:{parts[2]}"
                renamed_associations.append(association)
        self.go_associations = DataManager.create_annot_set_from_legacy_assocs(assocs=renamed_associations,
                                                                               ontology=self.go_ontology)

    def _load_expression_annotations(self, taxon_id: str, provider: str, source: str = "db"):
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

    def _load_do_annotations(self, taxon_id: str, provider: str, source: str = "db"):
        if source == "db":
            associations = []
            doa_list = get_disease_annotations(taxon_id=taxon_id)
            for doa in doa_list:
                associations.append(DataManager.create_annotation_record(
                    source_line="",
                    gene_id=doa["gene_id"],
                    gene_symbol=doa["gene_symbol"],
                    gene_type="gene",
                    taxon_id=taxon_id,
                    object_id=doa["do_id"],
                    qualifiers=[""],
                    aspect="D",
                    ecode="BMK" if doa["relationship_type"] == "is_marker_for" else "EXP",
                    references="",
                    prvdr=provider,
                    date=None))
            self.do_associations = DataManager.create_annot_set_from_legacy_assocs(
                assocs=associations, ontology=self.do_ontology)

    def load_annotations(self, associations_type: DataType, taxon_id: str, provider: str, source: str = "db"):
        if associations_type == DataType.GO:
            self._load_go_annotations(provider=provider)
        elif associations_type == DataType.EXPR:
            self._load_expression_annotations(taxon_id=taxon_id, provider=provider, source=source)
        elif associations_type == DataType.DO:
            self._load_do_annotations(taxon_id=taxon_id, source=source, provider=provider)
        return None

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

        node_type = ""
        if ontology_type == DataType.GO:
            node_type = "goterm"
        elif ontology_type == DataType.EXPR:
            node_type = "anatomyterm"
        if ontology_type == DataType.DO:
            node_type = "doterm"

        if source == "api":
            self._load_ontology_from_api(ontology, curie_prefix, node_type)
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
        elif ontology_type == DataType.DO:
            return "DOID"
        return ""

    def _load_ontology_from_api(self, ontology, curie_prefix: str, node_type: str):
        roots = get_ontology_roots(node_type=node_type)
        roots = [root for root in roots if root["curie"].startswith(curie_prefix)]
        visited_nodes = set(root["curie"] for root in roots)

        for root in roots:
            self.add_node_to_ontobio_ontology_if_not_exists(
                term_id=root["curie"],
                term_label=root["name"],
                term_type=root["namespace"],
                is_obsolete=False,
                ontology=ontology,
                check_exists=False)
        nodes = roots

        def process_node(node):
            new_children = []
            if node["descendantCount"] > 0:
                children = get_ontology_node_children(node_curie=node["curie"], node_type=node_type)
                for child in children:
                    if child["curie"] not in visited_nodes:
                        self.add_node_to_ontobio_ontology_if_not_exists(
                            term_id=child["curie"],
                            term_label=child["name"],
                            term_type=child["namespace"],
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
                    nodes.extend([node for node in future.result() if node.get("childCount", 0) > 0])

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
                                relation="subClassOf" if onto_pair["rel_type"].upper() == "IS_A" else "BFO:0000050")

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

    def get_best_human_orthologs(self, species_taxon: str, source: str = "db"):
        """Get best human orthologs for all genes from a given data provider."""
        if source not in ["db", "api"]:
            raise ValueError("source must be either 'db' or 'api'")
        if source == "api":
            raise NotImplementedError("API loading for human orthologs is not implemented yet")
        else:
            return get_best_human_orthologs_for_taxon(species_taxon)
