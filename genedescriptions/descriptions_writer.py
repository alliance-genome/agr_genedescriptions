import datetime
import json

from collections import OrderedDict
from typing import List

from genedescriptions.data_manager import DataManager
from genedescriptions.gene_description import GeneDescription
from genedescriptions.stats import DescriptionsOverallProperties, DescriptionsStats


class DescriptionsWriter(object):

    def __init__(self):
        self.overall_properties = DescriptionsOverallProperties()
        self.general_stats = DescriptionsStats()
        self.data = []

    def add_gene_desc(self, gene_description: GeneDescription):
        """add a gene description to the writer object

        Args:
            gene_description (GeneDescription): the gene description to be added
        """
        self.data.append(gene_description)

    def write_json(self, file_path: str, pretty: bool = False, include_single_gene_stats: bool = False,
                   data_manager: DataManager = None):
        """write the descriptions to a json file

        Args:
            file_path (str): the path to the file to write
            pretty (bool): whether to format the json file to make it more human-readable
            include_single_gene_stats (bool): whether to include statistics about the descriptions in the output file
            data_manager (DataManager): data manager containing the ontologies used to calculate onto stats
        """
        indent = 4 if pretty else None
        if include_single_gene_stats:
            for gene_desc in self.data:
                gene_desc.stats.calculate_stats(data_manager=data_manager)
            self.general_stats.calculate_stats(gene_descriptions=self.data)

        output = OrderedDict()
        output['overall_properties'] = vars(self.overall_properties)
        if include_single_gene_stats:
            output['general_stats'] = vars(self.general_stats)

        serialized_data = []
        for gd in self.data:
            entry = OrderedDict()
            entry['gene_id'] = gd.gene_id
            entry['gene_name'] = gd.gene_name
            entry['description'] = gd.description
            entry['go_description'] = gd.go_description
            entry['go_function_description'] = gd.go_function_description
            entry['go_process_description'] = gd.go_process_description
            entry['go_component_description'] = gd.go_component_description
            entry['do_description'] = gd.do_description
            entry['do_experimental_description'] = gd.do_experimental_description
            entry['do_biomarker_description'] = gd.do_biomarker_description
            entry['do_orthology_description'] = gd.do_orthology_description
            entry['orthology_description'] = gd.orthology_description
            entry['tissue_expression_description'] = gd.tissue_expression_description
            entry['gene_expression_cluster_description'] = gd.gene_expression_cluster_description
            entry['molecule_expression_cluster_description'] = gd.molecule_expression_cluster_description
            entry['anatomy_expression_cluster_description'] = gd.anatomy_expression_cluster_description
            entry['protein_domain_description'] = gd.protein_domain_description
            entry['human_gene_function_description'] = gd.human_gene_function_description
            entry['sister_species_description'] = gd.sister_species_description
            entry['publications'] = gd.publications
            entry['refs'] = gd.refs
            if include_single_gene_stats:
                entry['stats'] = vars(gd.stats)
            entry['paper_evidences'] = gd.paper_evidences
            entry['accession_evidences'] = gd.accession_evidences
            serialized_data.append(entry)
        output['data'] = serialized_data

        with open(file_path, "w") as outfile:
            json.dump(output, outfile, indent=indent)

        if include_single_gene_stats:
            for gene_desc in self.data:
                gene_desc.stats.delete_extra_info()

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
                                  release_version + " version of WormBase\"\n")
                    outfile.write("Automated_description\t\"" + genedesc.description +
                                  "\"\tPaper_evidence\t\"WBPaper00065943\"\n")
                    outfile.write("Automated_description\t\"" + genedesc.description +
                                  "\"\tPaper_evidence\t\"WBPaper00067038\"\n\n")

    def write_plain_text(self, file_path, header=None):
        """write the descriptions to a plain text file

        Args:
            file_path (str): the path to the file to write
            header (str): optional header to prepend to the file
        """
        with open(file_path, "w") as outfile:
            if header:
                outfile.write(header + "\n\n")
            for genedesc in self.data:
                if genedesc.description:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\n" + genedesc.description + "\n\n")
                else:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\nNo description available\n\n")

    def write_tsv(self, file_path, header=None):
        """write the descriptions to a tsv file

        Args:
            file_path (str): the path to the file to write
            header (str): optional header to prepend to the file
        """
        with open(file_path, "w") as outfile:
            if header:
                outfile.write(header + "\n\n")
            for genedesc in self.data:
                if genedesc.description:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\t" + genedesc.description + "\n")
                else:
                    outfile.write(genedesc.gene_id + "\t" + genedesc.gene_name + "\tNo description available\n")
