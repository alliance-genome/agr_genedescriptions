import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty

logger = logging.getLogger("Config Parser tests")


class TestConfigParser(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="INFO", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))

    def test_exclude_terms_list(self):
        self.assertTrue(all([term in ["GO:0008150", "GO:0003674", "GO:0005575", "GO:0005488", "GO:0005515",
                                      "GO:0044877"] for term in
                             self.conf_parser.get_module_property(module=Module.GO,
                                                                  prop=ConfigModuleProperty.EXCLUDE_TERMS)]),
                        "GO exclusion term list not loading")
        self.assertTrue(len(self.conf_parser.get_module_property(module=Module.DO_EXPERIMENTAL,
                                                                 prop=ConfigModuleProperty.EXCLUDE_TERMS)) > 0,
                        "DO terms exclusion not loading")

    def test_rename_terms(self):
        self.assertTrue(len(self.conf_parser.get_module_property(module=Module.GO,
                                                                 prop=ConfigModuleProperty.RENAME_TERMS)) == 7,
                        "GO term renaming list not loading")
        self.assertTrue(self.conf_parser.get_module_property(module=Module.DO_EXPERIMENTAL,
                                                             prop=ConfigModuleProperty.RENAME_TERMS) is None,
                        "DO term renaming list should be None")

    def test_evidence_codes(self):
        self.assertTrue("EXP" in list(self.conf_parser.get_evidence_codes_groups_map(module=Module.GO).keys()))
