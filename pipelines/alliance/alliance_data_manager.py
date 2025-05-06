import logging
from typing import List

from genedescriptions.commons import DataType, Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from pipelines.alliance.ateam_db_helper import get_expression_annotations, get_ontology_pairs, get_gene_data

logger = logging.getLogger(__name__)


class AllianceDataManager(DataManager):
    """
    This class is a subclass of DataManager and is used to manage data for the Alliance.
    It provides methods to load and save data, as well as to perform various data processing tasks.
    """

    def __init__(self, config: GenedescConfigParser, alliance_release_version: str = None):
        self.config = config
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

    def load_ontology_from_persistent_store(self, ontology_type: DataType, provider: str = None):
        curie_prefix = ""
        if ontology_type == DataType.GO:
            pass
        elif ontology_type == DataType.EXPR:
            if provider == "WB":
                curie_prefix = "WBbt"
            ontology = get_ontology_pairs(curie_prefix=curie_prefix)
            self.set_ontology(ontology_type=ontology_type, ontology=ontology, config=self.config)

    def load_gene_data_from_persistent_store(self, provider: str):
        self.gene_data = get_gene_data(provider=provider)
