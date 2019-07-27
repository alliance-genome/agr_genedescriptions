import datetime
import logging
import unittest
import os

from ontobio import AssociationSetFactory, OntologyFactory

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.ontology_tools import get_all_common_ancestors, find_set_covering, \
    set_all_information_content_values

logger = logging.getLogger("Gene Ontology Tools tests")


class TestOntologyTools(unittest.TestCase):

    def load_go_ontology(self):
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

    def load_do_ontology(self):
        logger.info("Starting Ontology Tools tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None)
        logger.info("Loading do ontology from file")
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        self.df.load_ontology_from_file(ontology_type=DataType.DO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "doid.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "doid.obo"),
                                        config=self.conf_parser)

    def test_get_common_ancestors(self):
        self.load_go_ontology()
        generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                              data_manager=self.df, config=self.conf_parser)
        node_ids = generator.terms_groups[('P', '')]["EXPERIMENTAL"]
        common_ancestors = get_all_common_ancestors(node_ids, generator.ontology)
        self.assertTrue(len(common_ancestors) > 0, "Common ancestors not found")
        associations = [association for subj_associations in self.df.go_associations.associations_by_subj.values() for
                        association in subj_associations]
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043055", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0061065", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043054", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043053", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        self.df.go_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                             ontology=self.df.go_ontology)
        self.conf_parser.config["go_sentences_options"]["exclude_terms"].append("GO:0040024")
        generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                              data_manager=self.df, config=self.conf_parser)
        node_ids = generator.terms_groups[('P', '')]["EXPERIMENTAL"]
        common_ancestors = get_all_common_ancestors(node_ids, generator.ontology)
        self.assertTrue("GO:0040024" not in common_ancestors, "Common ancestors contain blacklisted term")

    def test_information_content(self):

        #              0                   ic(0) = 0
        #            /| |\
        #           / | | \
        #          1  2 3  4               ic(1) = 0.693147181, ic(2) = 0.470003629, ic(3) = 0.980829253
        #         /\ /\/ \/
        #        /  5 6  7                 ic(5) = 0.980829253, ic(6) = 1.16315081, ic(7) = 1.16315081
        #       /  /\  \/
        #      /  8  9 10                  ic(8) = 1.049822124, ic(10) = 1.252762968
        #      \ / \/   \
        #      11  12   13                 ic(11) = 1.386294361, ic(12) = 1.386294361, ic(13) = 1.386294361

        ontology = OntologyFactory().create()
        for i in range(14):
            ontology.add_node(i, 'node' + str(i))
        ontology.add_parent(1, 0)
        ontology.add_parent(2, 0)
        ontology.add_parent(3, 0)
        ontology.add_parent(4, 0)
        ontology.add_parent(5, 1)
        ontology.add_parent(5, 2)
        ontology.add_parent(6, 2)
        ontology.add_parent(6, 3)
        ontology.add_parent(7, 3)
        ontology.add_parent(7, 4)
        ontology.add_parent(8, 5)
        ontology.add_parent(9, 5)
        ontology.add_parent(10, 6)
        ontology.add_parent(10, 7)
        ontology.add_parent(11, 1)
        ontology.add_parent(11, 8)
        ontology.add_parent(12, 8)
        ontology.add_parent(12, 9)
        ontology.add_parent(13, 10)
        set_all_information_content_values(ontology=ontology)
        self.assertTrue(ontology.node(0)["IC"] == 0, "Root IC not equal to 0")
        self.assertAlmostEqual(ontology.node(1)["IC"], 0.693147181)
        self.assertAlmostEqual(ontology.node(2)["IC"], 0.470003629)
        self.assertAlmostEqual(ontology.node(3)["IC"], 0.980829253)
        self.assertAlmostEqual(ontology.node(5)["IC"], 0.980829253)
        self.assertAlmostEqual(ontology.node(6)["IC"], 1.16315081)
        self.assertAlmostEqual(ontology.node(7)["IC"], 1.16315081)
        self.assertAlmostEqual(ontology.node(8)["IC"], 1.049822124)
        self.assertAlmostEqual(ontology.node(10)["IC"], 1.252762968)
        self.assertAlmostEqual(ontology.node(11)["IC"], 1.386294361)

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

    def test_set_covering_with_ontology(self):
        self.load_do_ontology()
        self.conf_parser.config["do_via_orth_sentences_options"]["trimming_algorithm"] = "ic"
        self.conf_parser.config["do_via_orth_sentences_options"]["max_num_terms"] = 5
        associations = [DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0080028", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0080056", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:14789", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0080026", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:14415", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0080045", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:3371", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:8886", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:674", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:5614", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:11830", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:8398", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:2256", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:5327", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:88452",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:1123", qualifiers="", aspect="D",
                                                             ecode="ISS", references="", prvdr="WB", date="")]
        self.df.do_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                             ontology=self.df.do_ontology)
        generator = OntologySentenceGenerator(gene_id="MGI:88452", module=Module.DO_ORTHOLOGY,
                                              data_manager=self.df, config=self.conf_parser)
        sentences = generator.get_module_sentences(
            aspect='D', qualifier='', merge_groups_with_same_prefix=True, keep_only_best_group=True,
            high_priority_term_ids=["DOID:0080028", "DOID:0080056", "DOID:14789", "DOID:0080026", "DOID:14415",
                                    "DOID:0080045"])
        print(sentences.get_description())

    def test_depth(self):
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "mgi_expr.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "mgi_expr.obo"),
                                        config=self.conf_parser)
        associations = [DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16748", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16469", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16668", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:18769", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32906", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16193", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16194", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16470", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32869", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17366", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17168", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32665", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:18215", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:18333", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16198", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17563", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32660", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32780", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16659", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17503", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17547", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17549", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32778", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16542", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17597", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17952", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16665", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:35164", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17550", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16728", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16727", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17536", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17458", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:18426", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17779", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32809", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17381", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17185", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17382", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:26720", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32845", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32767", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16092", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16894", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32842", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17787", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16896", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32763", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17544", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17021", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17540", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17606", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17382", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:35577", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17525", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16039", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:19157", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17577", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17412", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:18425", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16246", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32874", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32678", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32711", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32837", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16916", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16910", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16833", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17524", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:19155", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32672", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32870", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:32673", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:17373", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="MGI:101764",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="EMAPA:16089", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="MGI", date=""),
                        ]
        self.df.expression_associations = AssociationSetFactory().create_from_assocs(
            assocs=associations, ontology=self.df.expression_ontology)
        generator = OntologySentenceGenerator(gene_id="MGI:101764", module=Module.EXPRESSION,
                                              data_manager=self.df, config=self.conf_parser)
        sentences = generator.get_module_sentences(
            aspect='A', qualifier='Verified', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        self.assertTrue(sentences.get_description() == "is expressed in embryo")
