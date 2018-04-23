import unittest
import os

from config_parser import GenedescConfigParser
from data_fetcher import AGRDBDataFetcher
from descriptions_rules import generate_go_sentences
from neo4j.v1 import GraphDatabase

from descriptions_writer import GeneDesc


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config.yml"))
        self.exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.go_prepostfix_sentences_map = conf_parser.get_go_prepostfix_sentences_map()
        self.go_prepostfix_special_cases_sent_map = conf_parser.get_go_prepostfix_special_cases_sent_map()
        self.go_annotations_priority = conf_parser.get_go_annotations_priority()
        self.evidence_codes_groups_map = conf_parser.get_evidence_codes_groups_map()
        self.go_terms_exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.release_version = conf_parser.get_release("wb_data_fetcher")
        self.evidence_groups_priority_list = conf_parser.get_evidence_groups_priority_list()
        self.go_terms_replacement_dict = conf_parser.get_go_rename_terms()
        self.host = "localhost"
        self.port = "7687"
        self.uri = "bolt://" + self.host + ":" + self.port
        self.graph = GraphDatabase.driver(self.uri, auth=("neo4j", "neo4j"))
        self.exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.go_prepostfix_sentences_map = conf_parser.get_go_prepostfix_sentences_map()
        self.go_prepostfix_special_cases_sent_map = conf_parser.get_go_prepostfix_special_cases_sent_map()
        self.go_annotations_priority = conf_parser.get_go_annotations_priority()
        self.evidence_codes_groups_map = conf_parser.get_evidence_codes_groups_map()
        self.go_terms_exclusion_list = conf_parser.get_go_terms_exclusion_list()
        self.release_version = conf_parser.get_release("wb_data_fetcher")
        self.evidence_groups_priority_list = conf_parser.get_evidence_groups_priority_list()
        self.go_terms_replacement_dict = conf_parser.get_go_rename_terms()
        self.go_truncate_others_aggregation_word = conf_parser.get_go_truncate_others_aggregation_word()
        self.go_truncate_others_terms = conf_parser.get_go_truncate_others_terms()
        self.go_trim_min_distance_from_root = conf_parser.get_go_trim_min_distance_from_root()

    def test_read_from_neo4j(self):
        df = AGRDBDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              go_terms_replacement_dict=self.go_terms_replacement_dict, data_provider="WB",
                              db_graph=self.graph)
        for gene in df.get_gene_data():
            self.assertTrue(gene.id)

    def test_load_go_data(self):
        df = AGRDBDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              go_terms_replacement_dict=self.go_terms_replacement_dict, data_provider="WB",
                              db_graph=self.graph)
        df.load_go_data()
        ontology = df.get_go_ontology()
        # query alternative ID
        self.assertTrue(ontology.query_term("GO:0072661").id == "GO:0072659")
        parents = ["GO:1990778", "GO:0072657"]
        self.assertTrue(all(map(lambda x: x.id in parents, ontology.query_term("GO:0072661").get_parents())))

    def test_generate_descriptions(self):
        df = AGRDBDataFetcher(go_terms_exclusion_list=self.exclusion_list,
                              go_terms_replacement_dict=self.go_terms_replacement_dict, data_provider="WB",
                              db_graph=self.graph)
        df.load_go_data()
        counter = 0
        for gene in df.get_gene_data():
            gene_desc = GeneDesc(gene_id=gene.id, gene_name=gene.name)
            sentences = generate_go_sentences(df.get_go_annotations(
                gene.id, priority_list=self.go_annotations_priority, desc_stats=gene_desc.stats),
                go_ontology=df.get_go_ontology(),
                evidence_groups_priority_list=self.evidence_groups_priority_list,
                go_prepostfix_sentences_map=self.go_prepostfix_sentences_map,
                go_prepostfix_special_cases_sent_map=self.go_prepostfix_special_cases_sent_map,
                evidence_codes_groups_map=self.evidence_codes_groups_map,
                remove_parent_terms=True,
                merge_num_terms_threshold=3,
                merge_min_distance_from_root=self.go_trim_min_distance_from_root,
                desc_stats=gene_desc.stats, go_terms_replacement_dict=self.go_terms_replacement_dict,
                truncate_others_generic_word=self.go_truncate_others_aggregation_word,
                truncate_others_aspect_words=self.go_truncate_others_terms)
            if sentences:
                joined_sent = []
                func_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                    go_aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True,
                    desc_stats=gene_desc.stats)])
                if func_sent:
                    joined_sent.append(func_sent)
                proc_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                    go_aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True,
                    desc_stats=gene_desc.stats)])
                if proc_sent:
                    joined_sent.append(proc_sent)
                comp_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                    go_aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True,
                    desc_stats=gene_desc.stats)])
                if comp_sent:
                    joined_sent.append(comp_sent)

                go_desc = "; ".join(joined_sent) + "."
                if len(go_desc) > 0:
                    gene_desc.description = go_desc[0].upper() + go_desc[1:]
            else:
                gene_desc.description = "No description available"
            self.assertTrue(gene_desc.description)


