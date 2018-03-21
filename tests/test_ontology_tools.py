import unittest
import os

from config_parser import GenedescConfigParser
from data_fetcher import WBRawDataFetcher
from descriptions_rules import generate_go_sentences
from ontology_tools import calculate_nodes_distance_matrix, get_merged_term_ids_by_set_covering


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config.yml"))
        self.raw_files_source = conf_parser.get_raw_file_sources("wb_data_fetcher")
        self.chebi_file_source = conf_parser.get_chebi_file_source()
        self.cache_location = conf_parser.get_cache_location()
        self.species = conf_parser.get_wb_species()
        self.exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.go_prepostfix_sentences_map = conf_parser.get_go_prepostfix_sentences_map()
        self.go_prepostfix_special_cases_sent_map = conf_parser.get_go_prepostfix_special_cases_sent_map()
        self.go_annotations_priority = conf_parser.get_go_annotations_priority()
        self.evidence_codes_groups_map = conf_parser.get_evidence_codes_groups_map()
        self.go_terms_exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.release_version = conf_parser.get_release("wb_data_fetcher")
        self.evidence_groups_priority_list = conf_parser.get_evidence_groups_priority_list()
        self.go_terms_replacement_dict = conf_parser.get_go_rename_terms()

    def test_calculate_nodes_distance_matrix(self):
        df = WBRawDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              raw_files_source=self.raw_files_source,
                              chebi_file_url=self.chebi_file_source,
                              release_version=self.release_version,
                              species="c_elegans",
                              project_id=self.species["c_elegans"]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              go_terms_replacement_dict=self.go_terms_replacement_dict)
        df.load_go_data()
        sentences = []
        gene = next(df.get_gene_data())
        go_ids = [annotation["GO_ID"] for annotation in df.get_go_annotations(gene.id) if annotation["Aspect"] == 'P']
        dist_mat = calculate_nodes_distance_matrix(go_ids, df.get_go_ontology())
        self.assertTrue(True)

    def test_get_merged_term_ids_by_set_covering(self):
        df = WBRawDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              raw_files_source=self.raw_files_source,
                              chebi_file_url=self.chebi_file_source,
                              release_version=self.release_version,
                              species="c_elegans",
                              project_id=self.species["c_elegans"]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              go_terms_replacement_dict=self.go_terms_replacement_dict)
        df.load_go_data()
        get_merged_term_ids_by_set_covering()
        self.assertTrue(True)
