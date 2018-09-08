import datetime
import json
from typing import List

import numpy as np
import copy
from genedescriptions.descriptions_rules import DescriptionsStats, GeneDesc, DescriptionsOverallProperties


class DescriptionsWriter(object):

    def __init__(self):
        self.overall_properties = DescriptionsOverallProperties()
        self.general_stats = DescriptionsStats()
        self.data = []

    def add_gene_desc(self, gene_description: GeneDesc):
        """add a gene description to the writer object

        Args:
            gene_description (GeneDesc): the gene description to be added
        """
        self.data.append(gene_description)

    def _calculate_stats(self):
        """calculate overall stats and populate fields"""
        num_initial_go_ids_f_arr = [len(gene_desc.stats.set_initial_go_ids_f) for gene_desc in self.data if
                                    gene_desc.go_function_description is not None]
        self.general_stats.average_number_initial_go_terms_f = np.average(num_initial_go_ids_f_arr) if \
            len(num_initial_go_ids_f_arr) > 0 else 0
        num_initial_go_ids_p_arr = [len(gene_desc.stats.set_initial_go_ids_p) for gene_desc in self.data if
                                    gene_desc.go_process_description is not None]
        self.general_stats.average_number_initial_go_terms_p = np.average(num_initial_go_ids_p_arr) if \
            len(num_initial_go_ids_p_arr) > 0 else 0
        num_initial_go_ids_c_arr = [len(gene_desc.stats.set_initial_go_ids_c) for gene_desc in self.data if
                                    gene_desc.go_component_description is not None]
        self.general_stats.average_number_initial_go_terms_c = np.average(num_initial_go_ids_c_arr) if \
            len(num_initial_go_ids_c_arr) > 0 else 0
        num_final_go_ids_f_arr = [len(gene_desc.stats.set_final_go_ids_f) for gene_desc in self.data if
                                  gene_desc.go_function_description is not None]
        self.general_stats.average_number_final_go_terms_f = np.average(num_final_go_ids_f_arr) if \
            len(num_final_go_ids_f_arr) > 0 else 0
        num_final_go_ids_p_arr = [len(gene_desc.stats.set_final_go_ids_p) for gene_desc in self.data if
                                  gene_desc.go_process_description is not None]
        self.general_stats.average_number_final_go_terms_p = np.average(num_final_go_ids_p_arr) if \
            len(num_final_go_ids_p_arr) > 0 else 0
        num_final_go_ids_c_arr = [len(gene_desc.stats.set_final_go_ids_c) for gene_desc in self.data if
                                  gene_desc.go_component_description is not None]
        self.general_stats.average_number_final_go_terms_c = np.average(num_final_go_ids_c_arr) if \
            len(num_final_go_ids_c_arr) > 0 else 0
        num_initial_do_terms_arr = [len(gene_desc.stats.set_initial_do_ids) for gene_desc in self.data if
                                    gene_desc.do_description is not None]
        self.general_stats.average_number_initial_do_terms = np.average(num_initial_do_terms_arr) if \
            len(num_initial_do_terms_arr) > 0 else 0
        num_final_do_terms_arr = [len(gene_desc.stats.set_final_do_ids) for gene_desc in self.data if
                                  gene_desc.do_description is not None]
        self.general_stats.average_number_final_do_terms = np.average(num_final_do_terms_arr) if \
            len(num_final_do_terms_arr) > 0 else 0
        self.general_stats.total_number_of_genes = len(self.data)
        self.general_stats.number_genes_with_non_null_description = len([gene_desc for gene_desc in self.data if
                                                                         gene_desc.description is not None])
        self.general_stats.number_genes_with_non_null_go_description = len([gene_desc for gene_desc in self.data if
                                                                            gene_desc.go_description is not None])
        self.general_stats.number_genes_with_null_go_description = len([gene_desc for gene_desc in self.data if
                                                                        gene_desc.go_description is None])
        self.general_stats.number_genes_with_non_null_go_function_description = \
            len([gene_desc for gene_desc in self.data if gene_desc.go_function_description is not None])
        self.general_stats.number_genes_with_non_null_go_process_description = \
            len([gene_desc for gene_desc in self.data if gene_desc.go_process_description is not None])
        self.general_stats.number_genes_with_non_null_go_component_description = \
            len([gene_desc for gene_desc in self.data if gene_desc.go_component_description is not None])
        self.general_stats.number_genes_with_more_than_3_initial_go_terms = \
            len([gene_desc for gene_desc in self.data if len(gene_desc.stats.set_initial_go_ids_f) > 3 or
                 len(gene_desc.stats.set_initial_go_ids_p) > 3 or len(gene_desc.stats.set_initial_go_ids_c) > 3])
        self.general_stats.number_genes_with_non_null_do_description = len([gene_desc for gene_desc in self.data if
                                                                            gene_desc.do_description is not None])
        self.general_stats.number_genes_with_null_do_description = len([gene_desc for gene_desc in self.data if
                                                                        gene_desc.do_description is None])
        self.general_stats.number_genes_with_more_than_3_initial_do_terms = \
            len([gene_desc for gene_desc in self.data if len(gene_desc.stats.set_initial_do_ids) > 3])
        self.general_stats.number_genes_with_final_do_terms_covering_multiple_initial_terms = \
            len([gene_desc for gene_desc in self.data if
                 gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms > 0])
        num_go_annot_arr = [gene_desc.stats.total_number_go_annotations for gene_desc in self.data if
                            gene_desc.go_description is not None]
        self.general_stats.average_number_go_annotations = np.average(num_go_annot_arr) if len(num_go_annot_arr) > 0 \
            else 0
        num_do_annot_arr = [gene_desc.stats.total_number_do_annotations for gene_desc in self.data if
                            gene_desc.do_description is not None]
        self.general_stats.average_number_do_annotations = np.average(num_do_annot_arr) if len(num_do_annot_arr) > 0 \
            else 0
        self.general_stats.number_genes_with_more_than_3_best_orthologs = \
            len([gene_desc for gene_desc in self.data if len(gene_desc.stats.set_best_orthologs) > 3])
        self.general_stats.number_genes_with_non_null_orthology_description = len(
            [gene_desc for gene_desc in self.data if gene_desc.orthology_description is not None])
        self.general_stats.number_genes_with_null_orthology_description = len(
            [gene_desc for gene_desc in self.data if gene_desc.orthology_description is None])
        num_orthologs = [len(gene_desc.stats.set_best_orthologs) for gene_desc in self.data if
                         gene_desc.orthology_description is not None]
        self.general_stats.average_number_orthologs = np.average(num_orthologs) if len(num_orthologs) > 0 \
            else 0

    def write_json(self, file_path: str, pretty: bool = False, include_single_gene_stats: bool = False):
        """write the descriptions to a json file

        Args:
            file_path (str): the path to the file to write
            pretty (bool): whether to format the json file to make it more human-readable
            include_single_gene_stats (bool): whether to include statistics about the descriptions in the output file
        """
        indent = None
        if pretty:
            indent = 4
        if include_single_gene_stats:
            self._calculate_stats()
        json_serializable_self = copy.deepcopy(self)
        json_serializable_self.overall_properties = vars(json_serializable_self.overall_properties)
        if include_single_gene_stats:
            json_serializable_self.general_stats = vars(json_serializable_self.general_stats)
            for gene_desc in json_serializable_self.data:
                gene_desc.stats = vars(gene_desc.stats)
        else:
            for gene_desc in json_serializable_self.data:
                gene_desc.stats = None
            del json_serializable_self.general_stats
        json_serializable_self.data = [vars(gene_desc) for gene_desc in json_serializable_self.data]
        if not include_single_gene_stats:
            for gene_desc in json_serializable_self.data:
                del gene_desc["stats"]
        with open(file_path, "w") as outfile:
            json.dump(vars(json_serializable_self), outfile, indent=indent)

    def write_ace(self, file_path: str, curators_list: List[str], release_version: str):
        """write the descriptions to an ace file

        Args:
            file_path (str): the path to the file to write
            curators_list (List[str]): list of WBPerson Ids to be attached as evidences to the automated descriptions
        """
        now = datetime.datetime.now()
        with open(file_path, "w") as outfile:
            outfile.write("\n")
            for genedesc in self.data:
                if genedesc.description:
                    outfile.write("Gene : \"" + genedesc.gene_id[3:] + "\"\n")
                    # for evidence in genedesc.evidences:
                    #    accession_arr = evidence.split(":")
                    #    outfile.write("Automated_description\t\"" + genedesc.description +
                    #                  "\"\tAccession_evidence\t\"" + accession_arr[0] + "\" \"" + accession_arr[1] +
                    #                  "\"\n")
                    for curator in curators_list:
                        outfile.write("Automated_description\t\"" + genedesc.description +
                                      "\"\tCurator_confirmed\t\"" + curator + "\"\n")
                    # for paper in genedesc.papersref:
                    #    outfile.write("Automated_description\t\"" + genedesc.description +
                    #                  "\"\tPaper_evidence\t\"" + paper + "\"\n")
                    outfile.write("Automated_description\t\"" + genedesc.description + "\"\tDate_last_updated\t\"" +
                                  str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "\"\n")
                    outfile.write("Automated_description\t\"" + genedesc.description +
                                  "\"\tInferred_automatically\t\"" + "This description was generated automatically by a"
                                                                     " script based on data from the " +
                                  release_version + " version of WormBase\"\n\n")

    def write_plain_text(self, file_path):
        """write the descriptions to a plain text file

        Args:
            file_path (str): the path to the file to write
        """
        with open(file_path, "w") as outfile:
            for genedesc in self.data:
                if genedesc.description:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\n" + genedesc.description + "\n\n")
                else:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\nNo description available\n\n")

    def write_tsv(self, file_path):
        """write the descriptions to a tsv file

        Args:
            file_path (str): the path to the file to write
        """
        with open(file_path, "w") as outfile:
            for genedesc in self.data:
                if genedesc.description:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\t" + genedesc.description + "\n")
                else:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\tNo description available\n")
