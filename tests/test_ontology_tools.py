import logging
import math
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.ontology_tools import get_all_common_ancestors, find_set_covering, \
    set_all_information_content_values

logger = logging.getLogger("Gene Ontology Tools tests")


class TestOntologyTools(unittest.TestCase):

    def load_ontology(self):
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
            self.this_dir, "data", "gene_association_1.7.wb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.wb.partial"),
                                            config=self.conf_parser)

    def test_get_common_ancestors(self):
        self.load_ontology()
        generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                              data_manager=self.df, config=self.conf_parser)
        node_ids = generator.terms_groups[('P', '')]["EXPERIMENTAL"]
        common_ancestors = get_all_common_ancestors(node_ids, generator.ontology)
        self.assertTrue(len(common_ancestors) > 0, "Common ancestors not found")

    def test_information_content(self):
        self.load_ontology()
        set_all_information_content_values(ontology=self.df.go_ontology)
        roots = self.df.go_ontology.get_roots()
        for root_id in roots:
            self.assertTrue(self.df.go_ontology.node(root_id)["IC"] == 0, "Root IC not equal to 0")

    def test_find_set_covering(self):
        subsets = [("1", "1", {"A", "B", "C"}), ("2", "2", {"A", "B"}), ("3", "3", {"C"}), ("4", "4", {"A"}),
                   ("5", "5", {"B"}), ("6", "6", {"C"})]
        values = [2, 12, 5, 20, 20, 20]
        # test with weights
        set_covering = [best_set[0] for best_set in find_set_covering(subsets=subsets, value=values, max_num_subsets=3)]
        self.assertTrue("2" in set_covering)
        self.assertTrue("6" in set_covering)
        self.assertTrue("1" not in set_covering)
        self.assertTrue("3" not in set_covering)
        self.assertTrue("4" not in set_covering)
        self.assertTrue("5" not in set_covering)
        # test without weights
        set_covering_noweights = [best_set[0] for best_set in
                                  find_set_covering(subsets=subsets, value=None, max_num_subsets=3)]
        self.assertTrue("1" in set_covering_noweights and len(set_covering_noweights) == 1)
        # test wrong input
        costs_wrong = [1, 3]
        set_covering_wrong = find_set_covering(subsets=subsets, value=costs_wrong, max_num_subsets=3)
        self.assertTrue(set_covering_wrong is None, "Cost vector with length different than subsets should return None")

        subsets = [("1", "1", {"7"}), ("2", "2", {"7", "12", "13"}),
                   ("3", "3", {"16", "17"}), ("4", "4", {"11"}), ("6", "6", {"12", "13"}), ("7", "7", {"7"}),
                   ("9", "9", {"16", "17"}), ("11", "11", {"11"}), ("12", "12", {"12"}), ("13", "13", {"13"}),
                   ("16", "16", {"16"}), ("17", "17", {"17"})]
        values = [1, 1, 0.875061263, 1.301029996, 1.301029996, 1.602059991, 1.301029996, 1.698970004, 1.698970004,
                  1.698970004, 1.698970004, 1.698970004]
        set_covering = [best_set[0] for best_set in find_set_covering(subsets=subsets, value=values, max_num_subsets=3)]
        self.assertTrue(all([num in set_covering for num in ["2", "9", "11"]]))
