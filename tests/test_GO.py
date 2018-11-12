import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator

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
            self.this_dir, os.path.pardir, "tests", "data", "go.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, os.path.pardir, "tests",
                                                                         "cache", "go.obo"), config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, os.path.pardir, "tests", "data", "go_annotations.gaf"),
                                            associations_cache_path=os.path.join(self.this_dir, os.path.pardir, "tests",
                                                                                 "cache", "go_annotations.gaf"),
                                            config=self.conf_parser)

    def test_ontology_exists(self):
        self.assertTrue(self.df.go_ontology is not None)
        self.assertTrue(any(parent == "GO:0007052" for parent in
                            self.df.go_ontology.parents("GO:0000022")))

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

    def test_generate_sentence(self):
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000018", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("several processes" in sentences.get_description())
        self.assertTrue("DNA damage checkpoint" in sentences.get_description())
