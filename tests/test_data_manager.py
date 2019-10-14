import logging
import unittest
import os

from genedescriptions.commons import Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, DataType

logger = logging.getLogger("Gene Ontology Module tests")


class TestGOModule(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Starting DataManager tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading go ontology from file")
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "go_gd_test.obo"), ontology_cache_path=os.path.join(self.this_dir, "cache",
                                                                                       "go_gd_test.obo"),
                                        config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.wb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.wb.partial"),
                                            config=self.conf_parser)

    def test_ontology_exists(self):
        self.assertTrue(self.df.go_ontology is not None)
        self.assertTrue(any(parent == "GO:0009987" for parent in
                            self.df.go_ontology.parents("GO:0000075")))

    def test_annotations_exist(self):
        self.assertTrue(self.df.go_associations is not None)
        self.assertTrue(len(self.df.get_annotations_for_gene(
            gene_id="WB:WBGene00000001", annot_type=DataType.GO,
            include_obsolete=False, include_negative_results=False,
            priority_list=self.conf_parser.get_annotations_priority(module=Module.GO))) > 0)

    def test_rename_terms(self):
        self.assertTrue(all(len(self.df.go_ontology.search(term)) == 0 for term in list(
            self.conf_parser.get_module_property(module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS).keys())))

    def test_exclude_terms(self):
        test_annot = self.df.get_annotations_for_gene("WB:WBGene00000001", annot_type=DataType.GO)
        self.assertTrue(all([annot["object"]["id"] != "GO:0008286" for annot in test_annot]))

    def test_download_gz_file(self):
        test_file = self.df._get_cached_file(cache_path=os.path.join(self.this_dir, "cache",
                                                                     "c_elegans.PRJNA13758.WS273.geneIDs.txt.gz"),
                                             file_source_url="file://" + os.path.join(
                                                 self.this_dir, "data", "c_elegans.PRJNA13758.WS273.geneIDs.txt.gz"))
        self.assertTrue(test_file == os.path.join(self.this_dir, "cache", "c_elegans.PRJNA13758.WS273.geneIDs.txt"))

    def test_gene_data_functions(self):
        self.df.set_gene_data(gene_data=[Gene("1", "gene1", True, False), Gene("2", "gene2", False, True),
                                         Gene("3", "gene3", False, False), Gene("4", "gene4", True, True)])
        self.assertTrue(len([g for g in self.df.get_gene_data(include_dead_genes=False,
                                                              include_pseudo_genes=False)]) == 1)
        self.assertTrue(len([g for g in self.df.get_gene_data(include_dead_genes=True,
                                                              include_pseudo_genes=False)]) == 2)
        self.assertTrue(len([g for g in self.df.get_gene_data(include_dead_genes=False,
                                                              include_pseudo_genes=True)]) == 2)
        self.assertTrue(len([g for g in self.df.get_gene_data(include_dead_genes=True,
                                                              include_pseudo_genes=True)]) == 4)

    def test_get_human_gene_props(self):
        human_gene_props = self.df.get_human_gene_props()
        self.assertTrue(len(human_gene_props) > 0)

    def test_get_ensembl_hgnc_ids_map(self):
        ensembl_hgnc_ids_map = self.df.get_ensembl_hgnc_ids_map()
        self.assertTrue(len(ensembl_hgnc_ids_map) > 0)
