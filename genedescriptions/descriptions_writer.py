import datetime
import json
import copy

from collections import OrderedDict
from typing import List
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
            self.general_stats.calculate_stats(gene_descriptions=self.data)
            for gene_desc in self.data:
                gene_desc.stats.calculate_stats()
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
        for gene_desc in json_serializable_self.data:
            del gene_desc.add_gene_name
        json_serializable_self.data = [vars(gene_desc) for gene_desc in json_serializable_self.data]
        if not include_single_gene_stats:
            for gene_desc in json_serializable_self.data:
                del gene_desc["stats"]
        with open(file_path, "w") as outfile:
            json.dump(OrderedDict(vars(json_serializable_self)), outfile, indent=indent)

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
