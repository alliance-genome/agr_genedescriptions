import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBDataFetcher, DataType
from genedescriptions.descriptions_rules import generate_sentences


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config_wb.yml"))

    def test_generate_sentences(self):
        df = WBDataFetcher(raw_files_source=self.conf_parser.get_raw_file_sources("wb_data_fetcher"),
                           release_version="WS265", species="c_elegans",
                           project_id=self.conf_parser.get_wb_species()["c_elegans"]["project_id"],
                           cache_location=self.conf_parser.get_cache_location(), do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"])
        df.load_all_data_from_file(go_terms_replacement_regex=self.conf_parser.get_go_rename_terms(),
                                   go_terms_exclusion_list=self.conf_parser.get_go_terms_exclusion_list(),
                                   do_terms_replacement_regex=None,
                                   do_terms_exclusion_list=self.conf_parser.get_do_terms_exclusion_list())
        sentences = []
        for gene in df.get_gene_data():
            joined_sentence = []
            go_sent = generate_sentences(
                annotations=df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                        priority_list=self.conf_parser.get_go_annotations_priority()),
                evidence_groups_priority_list=self.conf_parser.get_go_evidence_groups_priority_list(),
                prepostfix_sentences_map=self.conf_parser.get_go_prepostfix_sentences_map(),
                prepostfix_special_cases_sent_map=self.conf_parser.get_go_prepostfix_special_cases_sent_map(),
                evidence_codes_groups_map=self.conf_parser.get_go_evidence_codes_groups_map(),
                ontology=df.go_ontology)
            if go_sent:
                joined_sentence.append("; ".join([sent.text for sent in go_sent.get_sentences(aspect='F')]))
            do_sent = generate_sentences(
                annotations=df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                        priority_list=self.conf_parser.get_do_annotations_priority()),
                evidence_groups_priority_list=self.conf_parser.get_do_evidence_groups_priority_list(),
                prepostfix_sentences_map=self.conf_parser.get_do_prepostfix_sentences_map(),
                evidence_codes_groups_map=self.conf_parser.get_do_evidence_codes_groups_map(),
                ontology=df.do_ontology)
            if do_sent:
                joined_sentence.append("; ".join([sent.text for sent in do_sent.get_sentences('D')]))
            sentences.append("; ".join(joined_sentence))
        self.assertGreater(len(sentences), 20000)

