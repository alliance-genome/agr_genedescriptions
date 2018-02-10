import json
import unittest
import os
from data_fetcher import WBRawDataFetcher
from descriptions_rules import generate_go_sentences


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        with open(os.path.join(this_dir, os.path.pardir, "config.json")) as conf_file:
            config = json.load(conf_file)
        self.raw_files_source = config["wb_data_fetcher_options"]["raw_files_source"]
        self.chebi_file_source = config["wb_data_fetcher_options"]["chebi_file_source"]
        self.cache_location = config["generic_data_fetcher_options"]["cache_location"]
        self.species = config["wb_data_fetcher_options"]["species"]
        self.evidence_codes = config["go_sentences_options"]["evidence_codes"]
        self.go_prepostfix_sentences_map = {(prepost["aspect"], prepost["group"]): (prepost["prefix"],
                                                                                    prepost["postfix"]) for
                                            prepost in config["go_sentences_options"]["go_prepostfix_sentences_map"]}
        self.go_annotations_priority = [name for name, priority in sorted([(ec["name"], ec["priority"]) for ec in
                                                                           self.evidence_codes], key=lambda x: x[1])]
        self.evidence_groups_list = list(set([evidence["group"] for evidence in self.evidence_codes]))
        self.evidence_codes_groups_map = {evidence["name"]: evidence["group"] for evidence in self.evidence_codes}

    def test_generate_go_sentences(self):
        df = WBRawDataFetcher(raw_files_source=self.raw_files_source, release_version="WS263",
                              species=self.species[3]["name"], project_id=self.species[3]["project_id"],
                              cache_location=self.cache_location, use_cache=False,
                              chebi_file_url=self.chebi_file_source)
        df.load_go_data()
        sentences = []
        for gene in df.get_gene_data():
            sentences.append(generate_go_sentences(df.get_go_annotations(gene.id), self.evidence_groups_list,
                                                   self.go_prepostfix_sentences_map, self.evidence_codes_groups_map))
        self.assertEqual(len(sentences), 48678)
