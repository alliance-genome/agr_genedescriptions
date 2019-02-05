import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, DataType

logger = logging.getLogger("Gene Ontology Module tests")


class TestGOModule(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading go ontology from file")
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "go_gd_test.obo"), ontology_cache_path=os.path.join(self.this_dir, "cache",
                                                                                       "go_gd_test.obo"),
                                        config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.wb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.wb.partial"),
                                            config=self.conf_parser)

    def test_ontology_exists(self):
        self.assertTrue(self.df.go_ontology is not None)
        self.assertTrue(any(parent == "GO:0009987" for parent in
                            self.df.go_ontology.parents("GO:0000075")))

    def test_annotations_exist(self):
        self.assertTrue(self.df.go_associations is not None)
        self.assertTrue(len(self.df.get_annotations_for_gene(
            gene_id="WB:WBGene00000001", annot_type=DataType.GO,
            include_obsolete=False, include_negative_results=False,
            priority_list=self.conf_parser.get_annotations_priority(module=Module.GO))) > 0)

    def test_rename_terms(self):
        self.assertTrue(all(len(self.df.go_ontology.search(term)) == 0 for term in list(
            self.conf_parser.get_module_property(module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS).keys())))

    def test_exclude_terms(self):
        pass
