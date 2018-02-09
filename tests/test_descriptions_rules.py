import configparser
import unittest
import os
from data_fetcher import WBRawDataFetcher, GO_ASPECT
from descriptions_rules import generate_go_sentences


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        config = configparser.ConfigParser()
        config.read(os.path.join(this_dir, os.path.pardir, "genedesc.ini"))
        self.raw_files_source = config.get("data_fetcher", "raw_files_source")
        self.cache_location = config.get("data_fetcher", "cache_location")
        self.species = config.get("generic", "species").split(",")[3]
        self.project_id = config.get("generic", "project_ids").split(",")[3]
        self.evidence_codes = config.get("go_sentences_options", "evidence_codes").split(",")
        self.evidence_codes = dict(zip(self.evidence_codes, range(len(self.evidence_codes))))
        self.evidence_groups = config.get("go_sentences_options", "evidence_groups").split(",")
        self.evidence_groups = dict(zip(self.evidence_groups, range(len(self.evidence_groups))))
        self.evidence_codes_groups_map = config.get("go_sentences_options", "evidence_codes_groups_map").split(",")
        self.evidence_codes_groups_map = dict(zip(self.evidence_codes, [int(num) for num in
                                                                        self.evidence_codes_groups_map]))
        self.go_prepostfix_sentences_list = config.get("go_sentences_options", "go_prepostfix_sentences_map").split(";")
        self.go_prepostfix_sentences_map = {}
        for go_prepostfix_sentence in self.go_prepostfix_sentences_list:
            indices, values = go_prepostfix_sentence.split(":")
            aspect, group = indices.split(",")
            prefix, postfix = values.split(",")
            self.go_prepostfix_sentences_map[(GO_ASPECT[aspect], self.evidence_groups[group])] = (prefix, postfix)

    def test_generate_go_sentences(self):
        df = WBRawDataFetcher(raw_files_source=self.raw_files_source, release_version="WS263", species=self.species,
                              project_id=self.project_id, cache_location=self.cache_location, use_cache=False)
        df.load_go_data()
        sentences = []
        for gene in df.get_gene_data():
            sentences.append(generate_go_sentences(df.get_go_annotations(gene.id), self.evidence_groups,
                                                   self.go_prepostfix_sentences_map, self.evidence_codes_groups_map))
        self.assertEqual(len(sentences), 48678)
