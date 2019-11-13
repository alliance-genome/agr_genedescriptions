import logging
import unittest
import os

from ontobio import AssociationSetFactory, OntologyFactory

from genedescriptions.commons import Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription
from genedescriptions.ontology_tools import set_ic_ontology_struct, get_all_common_ancestors, set_ic_annot_freq
from genedescriptions.precanned_modules import set_expression_module, set_gene_ontology_module
from wormbase.wb_data_manager import WBDataManager

logger = logging.getLogger(__name__)


class TestOntologyTools(unittest.TestCase):

    def setUp(self):
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "go_gd_test.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "go_gd_test.obo"),
                                        config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.wb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.wb.partial"),
                                            config=self.conf_parser)

    def test_trimming_lca(self):
        gene = Gene(id="WB:WBGene00000018", name="abl-1", dead=False, pseudo=False)
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "anatomy_ontology.WS274.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache",
                                                                         "anatomy_ontology.WS274.obo"),
                                        config=self.conf_parser)
        logger.info("Loading expression associations from file")
        self.conf_parser.config["expression_sentences_options"]["max_num_terms"] = 5
        self.conf_parser.config["expression_sentences_options"]["trim_min_distance_from_root"]["A"] = 4
        self.conf_parser.config["expression_sentences_options"]["remove_children_if_parent_is_present"] = False
        associations = [DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0006796", qualifiers=["Verified"],
                                                             aspect="A", ecode="IDA", references="", prvdr="WB",
                                                             date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0006759", qualifiers=["Verified"], aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0005300", qualifiers=["Verified"], aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0008598", qualifiers=["Verified"], aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0003681", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0005829", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0003927", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id=gene.id,
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="WBbt:0006751", qualifiers=["Verified"],
                                                             aspect="A",
                                                             ecode="IDA", references="", prvdr="WB", date=""),
                        ]
        self.df.expression_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                                     ontology=self.df.expression_ontology)
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "lca"
        self.conf_parser.config["expression_sentences_options"]["trimming_algorithm"] = "lca"
        gene_desc_lca = GeneDescription(gene_id=gene.id, config=self.conf_parser, gene_name="abl-1",
                                        add_gene_name=False)
        set_gene_ontology_module(dm=self.df, conf_parser=self.conf_parser, gene_desc=gene_desc_lca, gene=gene)
        set_expression_module(self.df, self.conf_parser, gene_desc_lca, gene)
        gene_desc_lca.stats.calculate_stats(data_manager=self.df)
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "ic"
        self.conf_parser.config["expression_sentences_options"]["trimming_algorithm"] = "ic"
        set_ic_ontology_struct(ontology=self.df.go_ontology, relations=self.df.go_relations)
        set_ic_ontology_struct(ontology=self.df.expression_ontology, relations=self.df.expr_relations)
        gene_desc_ic = GeneDescription(gene_id=gene.id, config=self.conf_parser, gene_name="abl-1",
                                       add_gene_name=False)
        set_gene_ontology_module(dm=self.df, conf_parser=self.conf_parser, gene_desc=gene_desc_ic, gene=gene)
        set_expression_module(self.df, self.conf_parser, gene_desc_ic, gene)
        gene_desc_ic.stats.calculate_stats(data_manager=self.df)
        self.assertTrue(gene_desc_lca.stats.coverage_percentage >= gene_desc_ic.stats.coverage_percentage)

