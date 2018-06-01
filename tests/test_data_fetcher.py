import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBRawDataFetcher, AGRRawDataFetcher


class TestRawDataFetcher(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config.yml"))
        raw_files_source = self.conf_parser.get_raw_file_sources("wb_data_fetcher")
        cache_location = self.conf_parser.get_cache_location()
        species = self.conf_parser.get_wb_species()
        self.exclusion_list = self.conf_parser.get_go_terms_exclusion_list()
        self.go_terms_replacement_dict = self.conf_parser.get_go_rename_terms()
        self.df = WBRawDataFetcher(raw_files_source=raw_files_source, release_version="WS265", species="c_elegans",
                                   project_id=species["c_elegans"]["project_id"],
                                   cache_location=cache_location, use_cache=False,
                                   go_terms_exclusion_list=self.exclusion_list,
                                   go_terms_replacement_dict=self.go_terms_replacement_dict)
        self.df.load_gene_data()
        self.df.load_go_data()
        self.df.load_disease_data()

    def test_get_parents(self):
        self.assertTrue(any(parent.id == "GO:0007052" for parent in
                            self.df.get_go_ontology().query_term("GO:0000022").get_parents()))

    def test_get_go_data(self):
        self.assertTrue(self.df.get_go_ontology())
        self.assertTrue(len(self.df.get_annotations_for_gene(geneid="WB:WBGene00000001")) > 0)

    def test_get_gene_data(self):
        self.assertGreater(len([gene for gene in self.df.get_gene_data()]), 20000)

    def test_data_fetchers(self):
        agr_df = AGRRawDataFetcher(raw_files_source=self.conf_parser.get_raw_file_sources("agr_data_fetcher"),
                                   cache_location=self.conf_parser.get_cache_location(), use_cache=False,
                                   go_terms_exclusion_list=self.exclusion_list,
                                   go_terms_replacement_dict=self.go_terms_replacement_dict,
                                   main_file_name=self.conf_parser.get_agr_mod_property("fb", "main_files"),
                                   bgi_file_name=self.conf_parser.get_agr_mod_property("fb", "bgi_file"),
                                   go_annotations_file_name=self.conf_parser.get_agr_mod_property("fb",
                                                                                                  "go_annotations"),
                                   organism_name=self.conf_parser.get_agr_mod_property("fb", "name"),
                                   release_version="1.0")
        agr_df.load_gene_data()
        agr_df.load_go_data()
        agr_df.load_disease_data()
        self.assertTrue(agr_df.get_go_ontology().query_term("GO:0003674"))


