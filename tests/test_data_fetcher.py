import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBRawDataFetcher


class TestWBRawDataFetcher(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config.yml"))
        raw_files_source = self.conf_parser.get_raw_file_sources("wb_data_fetcher")
        chebi_file_source = self.conf_parser.get_chebi_file_source()
        cache_location = self.conf_parser.get_cache_location()
        species = self.conf_parser.get_wb_species()
        exclusion_list = self.conf_parser.get_go_terms_exclusion_list()
        go_terms_replacement_dict = self.conf_parser.get_go_rename_terms()
        self.df = WBRawDataFetcher(raw_files_source=raw_files_source, release_version="WS263", species="c_elegans",
                                   project_id=species["c_elegans"]["project_id"],
                                   cache_location=cache_location, use_cache=False,
                                   chebi_file_url=chebi_file_source, go_terms_exclusion_list=exclusion_list,
                                   go_terms_replacement_dict=go_terms_replacement_dict)
        self.df.load_go_data()

    def test_get_parents(self):
        self.assertTrue(any(parent.id == "GO:0007052" for parent in
                            self.df.get_go_ontology().query_term("GO:0000022").get_parents()))

    def test_load_go_data(self):
        self.df.load_go_data()
        self.assertTrue(self.df.get_go_ontology())
        self.assertTrue(len(self.df.get_go_annotations(geneid="WBGene00000001")) > 0)
        

    def test_get_gene_data(self):
        self.df.load_go_data()
        self.assertEqual(len([gene for gene in self.df.get_gene_data()]), 48678)

