import unittest

from ontobio import OntologyFactory

from genedescriptions.commons import CommonAncestor
from genedescriptions.optimization import find_set_covering


class TestOptimization(unittest.TestCase):

    def test_find_set_covering(self):
        subsets = [CommonAncestor("1", "1", {"A", "B", "C"}), CommonAncestor("2", "2", {"A", "B"}),
                   CommonAncestor("3", "3", {"C"}), CommonAncestor("4", "4", {"A"}),
                   CommonAncestor("5", "5", {"B"}), CommonAncestor("6", "6", {"C"})]
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
        set_covering_noweights = [best_set[0] for best_set in find_set_covering(subsets=subsets, value=None,
                                                                                max_num_subsets=3)]
        self.assertTrue("1" in set_covering_noweights and len(set_covering_noweights) == 1)
        # test wrong input
        costs_wrong = [1, 3]
        set_covering_wrong = find_set_covering(subsets=subsets, value=costs_wrong, max_num_subsets=3)
        self.assertTrue(set_covering_wrong is None, "Cost vector with length different than subsets should return None")

        subsets = [CommonAncestor("1", "1", {"7"}), CommonAncestor("2", "2", {"7", "12", "13"}),
                   CommonAncestor("3", "3", {"16", "17"}), CommonAncestor("4", "4", {"11"}),
                   CommonAncestor("6", "6", {"12", "13"}), CommonAncestor("7", "7", {"7"}),
                   CommonAncestor("9", "9", {"16", "17"}), CommonAncestor("11", "11", {"11"}),
                   CommonAncestor("12", "12", {"12"}), CommonAncestor("13", "13", {"13"}),
                   CommonAncestor("16", "16", {"16"}), CommonAncestor("17", "17", {"17"})]
        values = [1, 1, 0.875061263, 1.301029996, 1.301029996, 1.602059991, 1.301029996, 1.698970004, 1.698970004,
                  1.698970004, 1.698970004, 1.698970004]
        set_covering = [best_set[0] for best_set in find_set_covering(subsets=subsets, value=values, max_num_subsets=3)]
        self.assertTrue(all([num in set_covering for num in ["2", "9", "11"]]))

    def test_set_covering_with_ontology(self):

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

        subsets = [CommonAncestor(node_id=1, node_label="1", covered_starting_nodes={"11", "12"}),
                   CommonAncestor(node_id=2, node_label="2", covered_starting_nodes={"11", "12", "13"}),
                   CommonAncestor(node_id=3, node_label="3", covered_starting_nodes={"13"}),
                   CommonAncestor(node_id=4, node_label="4", covered_starting_nodes={"13"}),
                   CommonAncestor(node_id=5, node_label="2", covered_starting_nodes={"11", "12"}),
                   CommonAncestor(node_id=6, node_label="6", covered_starting_nodes={"13"}),
                   CommonAncestor(node_id=7, node_label="7", covered_starting_nodes={"13"}),
                   CommonAncestor(node_id=8, node_label="8", covered_starting_nodes={"11", "12"}),
                   CommonAncestor(node_id=9, node_label="9", covered_starting_nodes={"12"}),
                   CommonAncestor(node_id=10, node_label="10", covered_starting_nodes={"13"}),
                   CommonAncestor(node_id=11, node_label="11", covered_starting_nodes={"11"}),
                   CommonAncestor(node_id=12, node_label="12", covered_starting_nodes={"12"}),
                   CommonAncestor(node_id=13, node_label="13", covered_starting_nodes={"13"})]

        values = [1, 1, 1, 1, 1, 1, 1, 20, 1, 1, 100, 1, 1]
        res = find_set_covering(subsets=subsets, ontology=ontology, value=values, max_num_subsets=2)
        self.assertTrue(all([sub[0] != 11 for sub in res]))
