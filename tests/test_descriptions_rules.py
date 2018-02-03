import configparser
import unittest
import os
from data_fetcher import WBRawDataFetcher
from descriptions_rules import generate_go_sentences, GOSentencesCollection


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        config = configparser.ConfigParser()
        config.read(os.path.join(this_dir, os.path.pardir, "genedesc.ini"))
        self.raw_files_source = config.get("data_fetcher", "raw_files_source")
        self.cache_location = config.get("data_fetcher", "cache_location")
        self.species = config.get("generic", "species").split(",")[3]
        self.project_id = config.get("generic", "project_ids").split(",")[3]

    def test_generate_go_sentences(self):
        df = WBRawDataFetcher(raw_files_source=self.raw_files_source, release_version="WS263", species=self.species,
                              project_id=self.project_id, cache_location=self.cache_location, use_cache=False)
        df.load_go_data()
        sentences = []
        for gene in df.get_gene_data():
            sentences.append(generate_go_sentences(df.get_go_annotations(gene.id)))
        self.assertEqual(len(sentences), 48678)
