import logging
import unittest
import os

from genedescriptions.commons import Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataType, DataManager
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription
from genedescriptions.precanned_modules import generate_ortholog_sentence_wormbase_human
from wormbase.wb_data_manager import WBDataManager

logger = logging.getLogger("Gene Ontology Module tests")


class TestGOModule(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, "config_test_wb.yml"))
        self.df = WBDataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"], config=self.conf_parser,
                                species="c_elegans")

    def test_load_expression_data(self):
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "anatomy_gd_test.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir,
                                                                         "cache", "anatomy_gd_test.obo"),
                                        config=self.conf_parser)
        self.df.load_associations_from_file(associations_type=DataType.EXPR, associations_url="file://" + os.path.join(
            self.this_dir, "data", "anatomy_gd_test.wb"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "anatomy_gd_test.wb"),
                                            config=self.conf_parser)
        self.assertTrue(self.df.expression_ontology is not None)
        self.assertTrue('WB:WBGene00000001' in self.df.expression_associations.associations_by_subj)
        for annotations in self.df.expression_associations.associations_by_subj.values():
            for annotation in annotations:
                self.assertTrue(annotation["evidence"]["type"] == "IDA")

    def test_load_disease_data(self):
        self.df.load_ontology_from_file(ontology_type=DataType.DO, ontology_url="file://" + os.path.join(
            self.this_dir, os.pardir, "data", "doid.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "doid.obo"),
                                        config=self.conf_parser)
        self.df.load_associations_from_file(associations_type=DataType.DO, associations_url=self.df.do_associations_url,
                                            associations_cache_path=os.path.join(self.this_dir, "cache", "do_ann.gaf"),
                                            association_additional_cache_path=os.path.join(self.this_dir, "cache",
                                                                                           "do_ann.daf"),
                                            association_additional_url=self.df.do_associations_new_url,
                                            config=self.conf_parser)
        self.assertTrue(any([annotation["evidence"]["type"] == "IMP" for annotations in
                             self.df.expression_associations.associations_by_subj.values() for annotation in
                             annotations]))

    def test_load_orthology_data(self):
        df = WBDataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"], config=self.conf_parser,
                           species="c_remanei")
        df.load_orthology_from_file()
        self.assertTrue(len(df.orthologs) > 0)

    def test_load_protein_domain_data(self):
        df = WBDataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"], config=self.conf_parser,
                           species="c_elegans")
        df.load_protein_domain_information()
        self.assertTrue(True)

    def test_expression_the_cell_renaming_to_widely(self):
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url=self.df.expression_ontology_url,
                                        ontology_cache_path=self.df.expression_ontology_cache_path,
                                        config=self.conf_parser)
        self.df.load_associations_from_file(associations_type=DataType.EXPR,
                                            associations_url=self.df.expression_associations_url,
                                            associations_cache_path=self.df.expression_associations_cache_path,
                                            config=self.conf_parser)
        gene_desc = GeneDescription(gene_id="WB:WBGene00007352", gene_name="cdc-48.1", add_gene_name=False)
        expr_sentence_generator = OntologySentenceGenerator(gene_id=gene_desc.gene_id, module=Module.EXPRESSION,
                                                            data_manager=self.df, config=self.conf_parser)
        expression_module_sentences = expr_sentence_generator.get_module_sentences(
            config=self.conf_parser, aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True,
            keep_only_best_group=False)
        gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=expression_module_sentences,
                                                                   module=Module.EXPRESSION)
        self.assertTrue("is expressed widely" in gene_desc.description)
