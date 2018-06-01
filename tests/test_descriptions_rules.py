import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBRawDataFetcher, AnnotationType
from genedescriptions.descriptions_rules import generate_sentences


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config_wb.yml"))
        self.raw_files_source = conf_parser.get_raw_file_sources("wb_data_fetcher")
        self.chebi_file_source = conf_parser.get_chebi_file_source()
        self.cache_location = conf_parser.get_cache_location()
        self.species = conf_parser.get_wb_species()
        self.exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.go_prepostfix_sentences_map = conf_parser.get_go_prepostfix_sentences_map()
        self.go_prepostfix_special_cases_sent_map = conf_parser.get_go_prepostfix_special_cases_sent_map()
        self.go_annotations_priority = conf_parser.get_go_annotations_priority()
        self.go_evidence_codes_groups_map = conf_parser.get_go_evidence_codes_groups_map()
        self.go_terms_exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.release_version = conf_parser.get_release("wb_data_fetcher")
        self.go_evidence_groups_priority_list = conf_parser.get_go_evidence_groups_priority_list()
        self.go_terms_replacement_dict = conf_parser.get_go_rename_terms()

        self.do_prepostfix_sentences_map = conf_parser.get_do_prepostfix_sentences_map()
        self.do_annotations_priority = conf_parser.get_do_annotations_priority()
        self.do_evidence_codes_groups_map = conf_parser.get_do_evidence_codes_groups_map()
        self.do_evidence_groups_priority_list = conf_parser.get_do_evidence_groups_priority_list()

    def test_generate_sentences(self):
        df = WBRawDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              raw_files_source=self.raw_files_source,
                              release_version=self.release_version,
                              species="c_elegans",
                              project_id=self.species["c_elegans"]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              go_terms_replacement_dict=self.go_terms_replacement_dict)
        df.load_go_data()
        df.load_disease_data()
        sentences = []
        for gene in df.get_gene_data():
            sentence = []
            go_sent = generate_sentences(annotations=df.get_annotations_for_gene(gene.id, annot_type=AnnotationType.GO),
                                         evidence_groups_priority_list=self.go_evidence_groups_priority_list,
                                         prepostfix_sentences_map=self.go_prepostfix_sentences_map,
                                         prepostfix_special_cases_sent_map=self.go_prepostfix_special_cases_sent_map,
                                         evidence_codes_groups_map=self.go_evidence_codes_groups_map,
                                         ontology=df.get_go_ontology())
            if go_sent:
                sentence.append(go_sent.get_sentences(aspect='F', ))
            do_sent = generate_sentences(annotations=df.get_annotations_for_gene(gene.id, annot_type=AnnotationType.DO),
                                         evidence_groups_priority_list=self.do_evidence_groups_priority_list,
                                         prepostfix_sentences_map=self.do_prepostfix_sentences_map,
                                         evidence_codes_groups_map=self.do_evidence_codes_groups_map,
                                         ontology=df.get_do_ontology())
            if do_sent:
                sentence.append(do_sent.get_sentences('D'))
        self.assertGreater(len(sentences), 20000)

