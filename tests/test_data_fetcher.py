import json
import unittest
import os
from data_fetcher import WBRawDataFetcher


class TestWBRawDataFetcher(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        with open(os.path.join(this_dir, os.path.pardir, "config.json")) as conf_file:
            config = json.load(conf_file)
        self.raw_files_source = config["wb_data_fetcher_options"]["raw_files_source"]
        self.chebi_file_source = config["wb_data_fetcher_options"]["chebi_file_source"]
        self.cache_location = config["generic_data_fetcher_options"]["cache_location"]
        self.species = config["wb_data_fetcher_options"]["species"]

    def test_get_gene_data(self):
        df = WBRawDataFetcher(raw_files_source=self.raw_files_source, release_version="WS263",
                              species=self.species[3]["name"], project_id=self.species[3]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              chebi_file_url=self.chebi_file_source)
        df.load_go_data()
        self.assertEqual(len([gene for gene in df.get_gene_data()]), 48678)

