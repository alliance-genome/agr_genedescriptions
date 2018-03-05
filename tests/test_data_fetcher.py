import unittest
import os

from config_parser import GenedescConfigParser
from data_fetcher import WBRawDataFetcher


class TestWBRawDataFetcher(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config.yml"))
        self.raw_files_source = conf_parser.get_raw_file_sources("wb_data_fetcher")
        self.chebi_file_source = conf_parser.get_chebi_file_source()
        self.cache_location = conf_parser.get_cache_location()
        self.species = conf_parser.get_wb_species()
        self.exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.go_terms_replacement_dict = conf_parser.get_go_rename_terms()

    def test_get_gene_data(self):
        df = WBRawDataFetcher(raw_files_source=self.raw_files_source, release_version="WS263",
                              species="c_elegans", project_id=self.species["c_elegans"]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              chebi_file_url=self.chebi_file_source, go_terms_exclusion_list=self.exclusion_list,
                              go_terms_replacement_dict=self.go_terms_replacement_dict)
        df.load_go_data()
        self.assertEqual(len([gene for gene in df.get_gene_data()]), 48678)

