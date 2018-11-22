import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from wormbase.wb_data_manager import WBDataManager

logger = logging.getLogger("Gene Ontology Module tests")


class TestGOModule(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, "config_test_wb.yml"))
        self.df = WBDataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"], config=self.conf_parser,
                                species="c_elegans")

    def test_load_expression_data(self):
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "anatomy_gd_test.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir,
                                                                         "cache", "anatomy_gd_test.obo"),
                                        config=self.conf_parser)
        self.df.load_associations_from_file(associations_type=DataType.EXPR, associations_url="file://" + os.path.join(
            self.this_dir, "data", "anatomy_gd_test.wb"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "anatomy_gd_test.wb"),
                                            config=self.conf_parser)
        self.assertTrue(self.df.expression_ontology is not None)
        self.assertTrue('WB:WBGene00000001' in self.df.expression_associations.associations_by_subj)
        for annotations in self.df.expression_associations.associations_by_subj.values():
            for annotation in annotations:
                self.assertTrue(annotation["evidence"]["type"] == "IDA")
