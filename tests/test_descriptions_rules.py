import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBRawDataFetcher
from genedescriptions.descriptions_rules import generate_go_sentences


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

    def test_generate_go_sentences(self):
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
        for gene in df.get_gene_data():
            sentences.append(generate_go_sentences(go_annotations=df.get_go_annotations(gene.id),
                                                   evidence_groups_priority_list=self.evidence_groups_priority_list,
                                                   go_prepostfix_sentences_map=self.go_prepostfix_sentences_map,
                                                   go_prepostfix_special_cases_sent_map=
                                                   self.go_prepostfix_special_cases_sent_map,
                                                   evidence_codes_groups_map=self.evidence_codes_groups_map,
                                                   go_ontology=df.get_go_ontology()))
        self.assertEqual(len(sentences), 48678)
