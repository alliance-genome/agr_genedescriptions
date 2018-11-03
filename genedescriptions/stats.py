import numpy as np


class SingleDescStats(object):
    """statistics for a single gene description"""

    def __init__(self):
        self.total_number_go_annotations = 0
        self.total_number_do_annotations = 0
        self.total_number_do_exp_bio_annotations = 0
        self.total_number_do_via_orth_annotations = 0
        self.number_final_do_term_covering_multiple_initial_do_terms = 0
        self.set_initial_experimental_go_ids_f = []
        self.set_initial_experimental_go_ids_p = []
        self.set_initial_experimental_go_ids_c = []
        self.set_initial_go_ids_f = []
        self.set_initial_go_ids_p = []
        self.set_initial_go_ids_c = []
        self.set_final_experimental_go_ids_f = []
        self.set_final_experimental_go_ids_p = []
        self.set_final_experimental_go_ids_c = []
        self.set_final_go_ids_f = []
        self.set_final_go_ids_p = []
        self.set_final_go_ids_c = []
        self.set_initial_do_ids = []
        self.set_final_do_ids = []
        self.set_best_orthologs = []
        self.set_initial_expression_ids = []
        self.set_final_expression_ids = []


class DescriptionsStats(object):
    """overall statistics for a set of gene descriptions"""

    def __init__(self):
        self.total_number_of_genes = 0
        self.number_genes_with_non_null_description = 0
        self.number_genes_with_non_null_go_description = 0
        self.number_genes_with_non_null_go_function_description = 0
        self.number_genes_with_non_null_go_process_description = 0
        self.number_genes_with_non_null_go_component_description = 0
        self.number_genes_with_null_go_description = 0
        self.number_genes_with_more_than_3_initial_go_terms = 0
        self.number_genes_with_non_null_do_description = 0
        self.number_genes_with_null_do_description = 0
        self.number_genes_with_more_than_3_initial_do_terms = 0
        self.number_genes_with_final_do_terms_covering_multiple_initial_terms = 0
        self.number_genes_with_non_null_do_experimental_description = 0
        self.number_genes_with_non_null_do_biomarker_description = 0
        self.number_genes_with_non_null_do_orthology_description = 0
        self.number_genes_with_non_null_protein_domain_description = 0
        self.number_genes_with_non_null_human_gene_function_description = 0
        self.number_genes_with_non_null_sister_species_description = 0
        self.number_genes_with_non_null_tissue_expression_description = 0
        self.number_genes_with_non_null_gene_expression_cluster_description = 0
        self.number_genes_with_non_null_molecule_expression_cluster_description = 0
        self.number_genes_with_non_null_anatomy_expression_cluster_description = 0
        self.average_number_initial_go_terms_f = 0
        self.average_number_initial_go_terms_p = 0
        self.average_number_initial_go_terms_c = 0
        self.average_number_final_go_terms_f = 0
        self.average_number_final_go_terms_p = 0
        self.average_number_final_go_terms_c = 0
        self.average_number_initial_do_terms = 0
        self.average_number_final_do_terms = 0
        self.average_number_go_annotations = 0
        self.average_number_do_annotations = 0
        self.average_number_orthologs = 0
        self.number_genes_with_non_null_orthology_description = 0
        self.number_genes_with_null_orthology_description = 0
        self.number_genes_with_more_than_3_best_orthologs = 0

    @staticmethod
    def _get_average_num_items_in_list_of_sets(set_var_name, desc_var_name, gene_descriptions):
        size_arr = [len(getattr(gene_desc.stats, set_var_name)) for gene_desc in gene_descriptions if
                    getattr(gene_desc, desc_var_name) is not None]
        return np.average(size_arr) if len(size_arr) > 0 else 0

    @staticmethod
    def _get_average(var_name, desc_var_name, gene_descriptions):
        clean_num_arr = [getattr(gene_desc.stats, var_name) for gene_desc in gene_descriptions if
                         getattr(gene_desc, desc_var_name)]
        return np.average(clean_num_arr) if len(clean_num_arr) > 0 else 0

    @staticmethod
    def _get_num_genes(gene_descriptions, desc_field_name: str, empty_desc: bool = False):
        return len([gene_desc for gene_desc in gene_descriptions if empty_desc and getattr(gene_desc, desc_field_name)
                    is None or not empty_desc and getattr(gene_desc, desc_field_name) is not None])

    def calculate_stats(self, gene_descriptions):
        """calculate overall stats and populate fields"""
        self.total_number_of_genes = len(gene_descriptions)
        self.average_number_initial_go_terms_f = self._get_average_num_items_in_list_of_sets(
            "set_initial_go_ids_f", "go_description", gene_descriptions)
        self.average_number_initial_go_terms_p = self._get_average_num_items_in_list_of_sets(
            "set_initial_go_ids_p", "go_description", gene_descriptions)
        self.average_number_initial_go_terms_c = self._get_average_num_items_in_list_of_sets(
            "set_initial_go_ids_c", "go_description", gene_descriptions)
        self.average_number_final_go_terms_f = self._get_average_num_items_in_list_of_sets(
            "set_final_go_ids_f", "go_description", gene_descriptions)
        self.average_number_final_go_terms_p = self._get_average_num_items_in_list_of_sets(
            "set_final_go_ids_p", "go_description", gene_descriptions)
        self.average_number_final_go_terms_c = self._get_average_num_items_in_list_of_sets(
            "set_final_go_ids_c", "go_description", gene_descriptions)
        self.average_number_initial_do_terms = self._get_average_num_items_in_list_of_sets(
            "set_initial_do_ids", "do_description", gene_descriptions)
        self.average_number_final_do_terms = self._get_average_num_items_in_list_of_sets(
            "set_final_do_ids", "do_description", gene_descriptions)
        self.number_genes_with_non_null_description = self._get_num_genes(gene_descriptions, "description")
        self.number_genes_with_non_null_go_description = self._get_num_genes(gene_descriptions, "go_description")
        self.number_genes_with_null_go_description = self._get_num_genes(gene_descriptions, "go_description", True)
        self.number_genes_with_non_null_go_function_description = self._get_num_genes(gene_descriptions,
                                                                                      "go_function_description")
        self.number_genes_with_non_null_go_process_description = self._get_num_genes(gene_descriptions,
                                                                                     "go_process_description")
        self.number_genes_with_non_null_go_component_description = self._get_num_genes(gene_descriptions,
                                                                                       "go_component_description")
        self.number_genes_with_more_than_3_initial_go_terms = \
            len([gene_desc for gene_desc in gene_descriptions if len(gene_desc.stats.set_initial_go_ids_f) > 3 or
                 len(gene_desc.stats.set_initial_go_ids_p) > 3 or len(gene_desc.stats.set_initial_go_ids_c) > 3])
        self.number_genes_with_non_null_do_description = self._get_num_genes(gene_descriptions, "do_description")
        self.number_genes_with_null_do_description = self._get_num_genes(gene_descriptions, "do_description", True)
        self.number_genes_with_non_null_do_experimental_description = self._get_num_genes(gene_descriptions,
                                                                                          "do_experimental_description")
        self.number_genes_with_non_null_do_biomarker_description = self._get_num_genes(gene_descriptions,
                                                                                       "do_biomarker_description")
        self.number_genes_with_non_null_do_orthology_description = self._get_num_genes(gene_descriptions,
                                                                                       "do_orthology_description")
        self.number_genes_with_non_null_tissue_expression_description = self._get_num_genes(
            gene_descriptions, "tissue_expression_description")
        self.number_genes_with_non_null_gene_expression_cluster_description = self._get_num_genes(
            gene_descriptions, "gene_expression_cluster_description")
        self.number_genes_with_non_null_molecule_expression_cluster_description = self._get_num_genes(
            gene_descriptions, "molecule_expression_cluster_description")
        self.number_genes_with_non_null_anatomy_expression_cluster_description = self._get_num_genes(
            gene_descriptions, "anatomy_expression_cluster_description")
        self.number_genes_with_non_null_protein_domain_description = self._get_num_genes(
            gene_descriptions, "protein_domain_description")
        self.number_genes_with_non_null_human_gene_function_description = self._get_num_genes(
            gene_descriptions, "human_gene_function_description")
        self.number_genes_with_non_null_sister_species_description = self._get_num_genes(gene_descriptions,
                                                                                         "sister_species_description")
        self.number_genes_with_more_than_3_initial_do_terms = \
            len([gene_desc for gene_desc in gene_descriptions if len(gene_desc.stats.set_initial_do_ids) > 3])
        self.number_genes_with_final_do_terms_covering_multiple_initial_terms = \
            len([gene_desc for gene_desc in gene_descriptions if
                 gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms > 0])
        self.average_number_go_annotations = self._get_average("total_number_go_annotations",
                                                               "go_description", gene_descriptions)
        self.average_number_do_annotations = self._get_average("total_number_do_annotations", "do_description",
                                                               gene_descriptions)
        self.number_genes_with_more_than_3_best_orthologs = \
            len([gene_desc for gene_desc in gene_descriptions if len(gene_desc.stats.set_best_orthologs) > 3])
        self.number_genes_with_non_null_orthology_description = self._get_num_genes(gene_descriptions,
                                                                                    "orthology_description")
        self.number_genes_with_null_orthology_description = self._get_num_genes(gene_descriptions,
                                                                                "orthology_description", True)
        self.average_number_orthologs = self._get_average_num_items_in_list_of_sets(
            "set_best_orthologs", "orthology_description", gene_descriptions)


class DescriptionsOverallProperties(object):
    def __init__(self, species: str = "", release_version: str = "", date: str = "", go_ontology_url: str = "",
                 go_association_url: str = "", do_ontology_url: str = "", do_association_url: str = ""):
        self.species = species
        self.release_version = release_version
        self.date = date
        self.go_ontology_url = go_ontology_url
        self.go_association_url = go_association_url
        self.do_ontology_url = do_ontology_url
        self.do_association_url = do_association_url

