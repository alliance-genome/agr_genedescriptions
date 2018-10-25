import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import WBDataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator


class TestDescriptionsRules(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config_wb.yml"))

    def test_generate_sentences(self):
        df = WBDataManager(raw_files_source=self.conf_parser.get_wb_raw_file_sources(),
                           release_version="WS265", species="c_elegans",
                           project_id=self.conf_parser.get_wb_organisms_info()["c_elegans"]["project_id"],
                           cache_location=self.conf_parser.get_cache_dir(),
                           do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"])
        df.load_all_data_from_file(go_terms_replacement_regex=self.conf_parser.get_module_simple_property(
            module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS),
                                   go_terms_exclusion_list=self.conf_parser.get_module_simple_property(
            module=Module.GO, prop=ConfigModuleProperty.EXCLUDE_TERMS),
                                   do_terms_replacement_regex=self.conf_parser.get_module_simple_property(
            module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.RENAME_TERMS),
                                   do_terms_exclusion_list=self.conf_parser.get_module_simple_property(
            module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.EXCLUDE_TERMS))
        sentences = []
        for gene in df.get_gene_data():
            go_sent_gen_common_props = self.conf_parser.get_sentence_generator_common_properties(Module.GO)
            go_sent_common_props = self.conf_parser.get_sentence_common_properties(module=Module.GO)
            annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                      priority_list=self.conf_parser.get_annotations_priority(
                                                          module=Module.GO))
            go_sent_generator = OntologySentenceGenerator(annotations=annotations, ontology=df.go_ontology,
                                                          **go_sent_gen_common_props)
            go_module_sentences = go_sent_generator.get_module_sentences(aspect='F', remove_parent_terms=True,
                                                                         **go_sent_common_props)
            if go_module_sentences.contains_sentences():
                sentences.append(go_module_sentences.get_description())
        self.assertTrue(any([True if sentence and len(sentence) > 0 else False for sentence in sentences]))

