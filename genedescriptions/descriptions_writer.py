import datetime
import json
from abc import ABCMeta, abstractmethod
import numpy as np
import copy
from genedescriptions.descriptions_rules import DescriptionsStats, GeneDesc


class DescriptionsWriter(metaclass=ABCMeta):

    @abstractmethod
    def __init__(self):
        self.general_stats = DescriptionsStats()
        self.data = []

    @abstractmethod
    def write(self):
        pass

    def add_gene_desc(self, gene_description: GeneDesc):
        """add a gene description to the writer object

        Args:
            gene_description (GeneDesc): the gene description to be added
        """
        self.data.append(gene_description)

    def _calculate_stats(self):
        """calculate overall stats and populate fields"""
        self.general_stats.average_number_initial_go_terms_f = np.average(
            [gene_desc.stats.number_initial_go_terms_f for gene_desc in self.data if gene_desc.go_function_description
             is not None])
        self.general_stats.average_number_initial_go_terms_p = np.average(
            [gene_desc.stats.number_initial_go_terms_p for gene_desc in self.data if gene_desc.go_process_description is
             not None])
        self.general_stats.average_number_initial_go_terms_c = np.average(
            [gene_desc.stats.number_initial_go_terms_c for gene_desc in self.data if gene_desc.go_component_description
             is not None])
        self.general_stats.average_number_final_go_terms_f = np.average(
            [gene_desc.stats.number_final_go_terms_f for gene_desc in self.data if gene_desc.go_function_description is
             not None])
        self.general_stats.average_number_final_go_terms_p = np.average(
            [gene_desc.stats.number_final_go_terms_p for gene_desc in self.data if gene_desc.go_process_description is
             not None])
        self.general_stats.average_number_final_go_terms_c = np.average(
            [gene_desc.stats.number_final_go_terms_c for gene_desc in self.data if gene_desc.go_component_description is
             not None])
        self.general_stats.average_number_initial_do_terms = np.average(
            [gene_desc.stats.number_initial_do_terms for gene_desc in self.data if gene_desc.do_description is not
             None])
        self.general_stats.average_number_final_do_terms = np.average(
            [gene_desc.stats.number_final_do_terms for gene_desc in self.data if gene_desc.do_description is not None])
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
            len([gene_desc for gene_desc in self.data if gene_desc.stats.number_initial_go_terms_f > 3 or
                 gene_desc.stats.number_initial_go_terms_p > 3 or gene_desc.stats.number_initial_go_terms_c > 3])
        self.general_stats.number_genes_with_non_null_do_description = len([gene_desc for gene_desc in self.data if
                                                                            gene_desc.do_description is not None])
        self.general_stats.number_genes_with_null_do_description = len([gene_desc for gene_desc in self.data if
                                                                        gene_desc.do_description is None])
        self.general_stats.number_genes_with_more_than_3_initial_do_terms = \
            len([gene_desc for gene_desc in self.data if gene_desc.stats.number_initial_do_terms > 3])
        self.general_stats.number_genes_with_final_do_terms_covering_multiple_initial_terms = \
            len([gene_desc for gene_desc in self.data if
                 gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms > 0])
        self.general_stats.average_number_go_annotations = np.average(
            [gene_desc.stats.total_number_go_annotations for gene_desc in self.data if gene_desc.go_description is not
             None])
        self.general_stats.average_number_do_annotations = np.average(
            [gene_desc.stats.total_number_do_annotations for gene_desc in self.data if gene_desc.do_description is not
             None])


class JsonGDWriter(DescriptionsWriter):
    """generate gene descriptions in json format"""
    def __init__(self):
        super().__init__()

    def write(self, file_path: str, pretty: bool = False, include_single_gene_stats: bool = False):
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


class WBWriter(DescriptionsWriter):
    def __init__(self):
        super().__init__()

    def write(self, file_path: str):
        """write the descriptions to a WB file

        Args:
            file_path (str): the path to the file to write
        """
        with open(file_path, "w") as outfile:
            for genedesc in self.data:
                now = datetime.datetime.now()
                outfile.write(genedesc.gene_id + "\t" + str(now.year) + "-" + str(now.month) + "-" + str(now.day) +
                              "\t" + genedesc.publications + "\t" + genedesc.refs + "\t" + genedesc.description + "\t" +
                              genedesc.species + "\t" + "This description was generated automatically by a script "
                                                        "based on homology/orthology data, Gene Ontology (GO) "
                                                        "annotations, Disease ontology (DO) annotations, and tissue "
                                                        "expression data from the " + genedesc.release_version +
                              " version of WormBase)")
