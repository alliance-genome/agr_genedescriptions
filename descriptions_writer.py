import json
from abc import ABCMeta, abstractmethod
from collections import defaultdict
from typing import List
import numpy as np
import copy

from data_fetcher import AGRDBDataFetcher


class SingleDescStats(object):
    """statistics for a single gene description"""
    def __init__(self):
        self.num_terms_notrim_nogroup_priority_nomerge = defaultdict(int)
        self.num_terms_trim_nogroup_priority_nomerge = defaultdict(int)
        self.num_terms_trim_group_priority_merge = defaultdict(int)
        self.total_num_go_annotations = 0
        self.num_prioritized_go_annotations = 0
        self.terms_notrim_nogroup_priority_nomerge = defaultdict(list)
        self.terms_trim_nogroup_priority_nomerge = defaultdict(list)
        self.terms_trim_group_priority_merge = defaultdict(list)


class GeneDesc(object):
    """gene description"""
    def __init__(self, gene_id: str, gene_name: str = "", description: str = "", stats: SingleDescStats = None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = description
        if stats:
            self.stats = stats
        else:
            self.stats = SingleDescStats()


class DescriptionsStats(object):
    """overall statistics for a set of gene descriptions"""
    def __init__(self):
        self.num_genes_with_go_sentence = 0
        self.average_num_go_terms_if_desc_trim_group_priority_merge = 0
        self.average_num_go_terms_if_desc_trim_nogroup_priority_nomerge = 0
        self.average_num_go_terms_if_desc_notrim_nogroup_priority_nomerge = 0


class DescriptionsWriter(metaclass=ABCMeta):

    @abstractmethod
    def __init__(self):
        self.data = []
        self.general_stats = DescriptionsStats()

    @abstractmethod
    def write(self):
        pass

    def add_gene_desc(self, gene_description: GeneDesc):
        """add a gene description to the writer object

        :param gene_description: the gene description to be added
        :type gene_description: GeneDesc
        """
        self.data.append(gene_description)

    def _calculate_stats(self):
        """calculate overall stats and populate fields"""
        self.general_stats.average_num_go_terms_if_desc_trim_group_priority_merge = np.average(
            [sum(gene_desc.stats.num_terms_trim_group_priority_merge.values()) for
             gene_desc in self.data if gene_desc.description != "No description available"])
        self.general_stats.average_num_go_terms_if_desc_trim_nogroup_priority_nomerge = np.average(
            [sum(gene_desc.stats.num_terms_trim_nogroup_priority_nomerge.values()) for gene_desc in self.data if
             gene_desc.description != "No description available"])
        self.general_stats.average_num_go_terms_if_desc_notrim_nogroup_priority_nomerge = np.average(
            [sum(gene_desc.stats.num_terms_notrim_nogroup_priority_nomerge.values()) for gene_desc in self.data if
             gene_desc.description != "No description available"])
        self.general_stats.num_genes_with_go_sentence = len([gene_desc for gene_desc in self.data if
                                                             gene_desc.description != "No description available"])


class JsonGDWriter(DescriptionsWriter):
    """generate gene descriptions in json format"""
    def __init__(self):
        super().__init__()

    def write(self, file_path: str, pretty: bool = False, include_single_gene_stats: bool = False):
        """write the descriptions to a json file

        :param file_path: the path to the file to write
        :type file_path: str
        :param pretty: whether to format the json file to make it more human-readable
        :type pretty: bool
        :param include_single_gene_stats: whether to include statistics about the descriptions in the output file
        :type include_single_gene_stats: bool
        """
        indent = None
        if pretty:
            indent = 4
        self._calculate_stats()
        json_serializable_self = copy.deepcopy(self)
        json_serializable_self.general_stats = vars(json_serializable_self.general_stats)
        if include_single_gene_stats:
            for gene_desc in json_serializable_self.data:
                gene_desc.stats = vars(gene_desc.stats)
        else:
            for gene_desc in json_serializable_self.data:
                gene_desc.stats = None
        json_serializable_self.data = [vars(gene_desc) for gene_desc in json_serializable_self.data]
        if not include_single_gene_stats:
            for gene_desc in json_serializable_self.data:
                del gene_desc["stats"]
        with open(file_path, "w") as outfile:
            json.dump(vars(json_serializable_self), outfile, indent=indent)


class Neo4jGDWriter(DescriptionsWriter):
    """write gene descriptions to AGR neo4j database"""

    def __init__(self):
        super().__init__()

    def write(self, db_graph):
        query = """
            UNWIND $descriptions as row 

            MATCH (g:GOTerm:Ontology {primaryKey:row.gene_id})
                SET g.automatedGeneSynopsis = row.description
            """
        AGRDBDataFetcher.query_db(db_graph=db_graph, query=query, parameters={"descriptions": self.data})

