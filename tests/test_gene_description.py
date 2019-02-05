import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription

logger = logging.getLogger("Gene Ontology Module tests")


class TestDescriptionsGenerator(unittest.TestCase):

    def setUp(self):
        logger.info("Starting Ontology Tools tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading go ontology from file")
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "go_gd_test.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "go_gd_test.obo"),
                                        config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.fb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.fb.partial"),
                                            config=self.conf_parser)
        logging.basicConfig(filename=None, level="INFO", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')

    def test_set_or_extend_module_description_and_final_stats(self):
        gene_desc = GeneDescription(gene_id="FB:FBgn0027655", gene_name="Test gene", add_gene_name=False)
        go_sent_generator = OntologySentenceGenerator(gene_id="FB:FBgn0027655", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.GO_PROCESS, module_sentences=sentences)
        self.assertTrue(gene_desc.description, "Is involved in several processes, including axo-dendritic transport, "
                                               "establishment of mitotic spindle orientation, and positive regulation "
                                               "of extent of heterochromatin assembly")
        gene_desc = GeneDescription(gene_id="FB:FBgn0027655", gene_name="Test gene", add_gene_name=True)
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.GO_PROCESS, module_sentences=sentences)
        self.assertTrue(gene_desc.description, "Test gene is involved in several processes, including axo-dendritic "
                                               "transport, establishment of mitotic spindle orientation, and positive "
                                               "regulation of extent of heterochromatin assembly")


