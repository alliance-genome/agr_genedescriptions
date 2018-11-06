import logging
import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager, DataType

logger = logging.getLogger("Gene Descriptions tests")


class TestRawDataFetcher(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])

    def test_load_from_file(self):
        logger.info("Testing loading go ontology from file")
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, os.path.pardir, "tests", "data", "go.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, os.path.pardir, "tests",
                                                                         "cache", "go.obo"), config=self.conf_parser)
        self.assertTrue(self.df.go_ontology is not None)
        self.assertTrue(any(parent == "GO:0007052" for parent in
                            self.df.go_ontology.parents("GO:0000022")))
        logger.info("Testing loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, os.path.pardir, "tests", "data", "go_annotations.gaf"),
                                            associations_cache_path=os.path.join(self.this_dir, os.path.pardir, "tests",
                                                                                 "cache", "go_annotations.gaf"),
                                            config=self.conf_parser)
        self.assertTrue(self.df.go_associations is not None)
