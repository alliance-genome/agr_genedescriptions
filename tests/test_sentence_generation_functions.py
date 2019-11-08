import os
import unittest

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.sentence_generation_functions import rename_human_ortholog_name, compose_sentence


class TestSentenceGenerationFunctions(unittest.TestCase):

    def test_rename_human_ortholog_name(self):
        self.assertTrue("family member" not in rename_human_ortholog_name("nuclear pore complex interacting protein "
                                                                          "family member A8"))

    def test_compose_sentence(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "tests", "config_test.yml"))
        sentence = compose_sentence(prefix="Is expressed in", additional_prefix="several processes, including",
                                    term_names=["cell", "tail", "head", "male"],
                                    postfix="based on experimental observation",
                                    config=conf_parser, ancestors_with_multiple_children={"head"}, rename_cell=True,
                                    put_anatomy_male_at_end=True)
        self.assertTrue("cell" not in sentence)
        self.assertTrue("and in male" in sentence)
        sentence = compose_sentence(prefix="Is expressed in", additional_prefix="several processes, including",
                                    term_names=["cell"],
                                    postfix="based on experimental observation",
                                    config=conf_parser, rename_cell=True)
        self.assertTrue(sentence == "Is expressed widely based on experimental observation")

