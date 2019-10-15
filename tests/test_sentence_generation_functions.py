import unittest

from genedescriptions.sentence_generation_functions import rename_human_ortholog_name


class TestSentenceGenerationFunctions(unittest.TestCase):

    def test_rename_human_ortholog_name(self):
        self.assertTrue("family member" not in rename_human_ortholog_name("nuclear pore complex interacting protein "
                                                                          "family member A8"))
