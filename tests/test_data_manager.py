import unittest
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import WBDataManager, DataType


class TestRawDataFetcher(unittest.TestCase):

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "config_wb.yml"))
        species = self.conf_parser.get_wb_species()
        self.df = WBDataManager(raw_files_source=self.conf_parser.get_raw_file_sources("wb_data_fetcher"),
                                release_version="WS266", species="c_elegans",
                                project_id=species["c_elegans"]["project_id"],
                                cache_location=self.conf_parser.get_cache_location(), do_relations=None,
                                go_relations=["subClassOf", "BFO:0000050"])

    def test_load_gene_data_from_file(self):
        self.df.load_gene_data_from_file()
        self.assertGreater(len([gene for gene in self.df.get_gene_data()]), 20000)

    def test_load_go_ontology_from_file(self):
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=self.df.go_ontology_url,
                                        ontology_cache_path=self.df.go_ontology_cache_path,
                                        terms_replacement_regex=self.conf_parser.get_go_rename_terms())
        self.assertTrue(self.df.go_ontology is not None)
        self.assertTrue(any(parent == "GO:0007052" for parent in
                            self.df.go_ontology.parents("GO:0000022")))

    def test_load_go_associations_from_file(self):
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=self.df.go_ontology_url,
                                        ontology_cache_path=self.df.go_ontology_cache_path,
                                        terms_replacement_regex=self.conf_parser.get_go_rename_terms())
        self.df.load_associations_from_file(associations_type=DataType.GO,
                                            associations_url=self.df.go_associations_url,
                                            associations_cache_path=self.df.go_associations_cache_path,
                                            exclusion_list=self.conf_parser.get_go_terms_exclusion_list())
        self.assertTrue(self.df.go_associations is not None)
        self.assertTrue(len(self.df.get_annotations_for_gene(gene_id="WB:WBGene00000001")) > 0)

    def test_load_do_ontology_from_file(self):
        self.df.load_ontology_from_file(ontology_type=DataType.DO, ontology_url=self.df.do_ontology_url,
                                        ontology_cache_path=self.df.do_ontology_cache_path)
        self.assertTrue(self.df.do_ontology is not None)

    def test_load_do_associations_from_file(self):
        self.df.load_ontology_from_file(ontology_type=DataType.DO, ontology_url=self.df.do_ontology_url,
                                        ontology_cache_path=self.df.do_ontology_cache_path)
        self.df.load_associations_from_file(associations_type=DataType.DO,
                                            associations_url=self.df.do_associations_url,
                                            associations_cache_path=self.df.do_associations_cache_path,
                                            exclusion_list=self.conf_parser.get_do_terms_exclusion_list())
        self.assertTrue(self.df.do_associations is not None)

    def test_load_expression_ontology_from_file(self):
        self.df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url=self.df.expression_ontology_url,
                                        ontology_cache_path=self.df.expression_ontology_cache_path,
                                        terms_replacement_regex=self.conf_parser.get_expression_rename_terms())
        self.assertTrue(self.df.expression_ontology is not None)

    def test_load_orthology_from_file(self):
        species = self.conf_parser.get_wb_species()
        df = WBDataManager(raw_files_source=self.conf_parser.get_raw_file_sources("wb_data_fetcher"),
                           release_version="WS265", species="c_briggsae",
                           project_id=species["c_briggsae"]["project_id"],
                           cache_location=self.conf_parser.get_cache_location(), do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"], sister_sp_fullname="Caenorhabditis elegans")
        sister_df = WBDataManager(raw_files_source=self.conf_parser.get_raw_file_sources("wb_data_fetcher"),
                                  release_version="WS265", species="c_elegans",
                                  project_id=species["c_elegans"]["project_id"],
                                  cache_location=self.conf_parser.get_cache_location(), do_relations=None,
                                  go_relations=["subClassOf", "BFO:0000050"])
        sister_df.load_gene_data_from_file()
        sister_df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=sister_df.go_ontology_url,
                                          ontology_cache_path=sister_df.go_ontology_cache_path,
                                          terms_replacement_regex=self.conf_parser.get_go_rename_terms())
        sister_df.load_associations_from_file(associations_type=DataType.GO,
                                              associations_url=sister_df.go_associations_url,
                                              associations_cache_path=sister_df.go_associations_cache_path,
                                              exclusion_list=self.conf_parser.get_do_terms_exclusion_list())
        df.load_orthology_from_file()
        best_orthologs, curr_orth_fullname = df.get_best_orthologs_for_gene(
            "WB:WBGene00000307", ["Caenorhabditis elegans"], sister_df, ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                                                         "HTP", "HDA", "HMP", "HGI", "HEP"])
        self.assertTrue(best_orthologs)


