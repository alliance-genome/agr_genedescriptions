import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from wormbase.wb_data_manager import WBDataManager


logger = logging.getLogger("Gene Descriptions tests")


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="DEBUG", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting Gene Description generation tests")
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "wormbase", "config_wb.yml"))

    def test_generate_sentences(self):
        logger.info("Testing generating descriptions")
        dm = WBDataManager(raw_files_source=self.conf_parser.get_wb_raw_file_sources(),
                           release_version="WS265", species="c_elegans",
                           project_id=self.conf_parser.get_wb_organisms_info()["c_elegans"]["project_id"],
                           cache_location=self.conf_parser.get_cache_dir(),
                           do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"])
        dm.load_all_data_from_file(config=self.conf_parser)
        sentences = []
        for gene in dm.get_gene_data():
            go_sent_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.GO, data_manager=dm,
                                                          config=self.conf_parser)
            go_module_sentences = go_sent_generator.get_module_sentences(aspect='F', config=self.conf_parser)
            if go_module_sentences.contains_sentences():
                sentences.append(go_module_sentences.get_description())
        self.assertTrue(any([True if sentence and len(sentence) > 0 else False for sentence in sentences]))

