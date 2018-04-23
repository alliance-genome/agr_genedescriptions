import gzip
import json
import tarfile
import urllib.request
import shutil
import os
import logging
import re
from itertools import chain

import goatools
from goatools.obo_parser import GODag
from Bio.UniProt.GOA import gafiterator
from abc import ABCMeta, abstractmethod
from collections import namedtuple, defaultdict
from typing import List, Iterable, Dict, Tuple

from descriptions_writer import SingleDescStats

Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


def get_parents(self):
    """Return parent GO IDs."""
    return set([parent for parent in chain(self.parents, getattr(self, "relationship", defaultdict(set))["part_of"])])


class GOTerm(object):
    """go term with the same properties and methods defined by goatools GOTerm

    only the properties used for gene descriptions are implemented
    """

    def __init__(self, name: str, depth: int, node_id: str, parents: List, children: List, ontology):
        self.depth = depth
        self.name = name
        self.id = node_id
        self._parents = parents
        self._children = children
        self._ontology = ontology

    def get_parents(self) -> List["GOTerm"]:
        """get the parent terms of the current term

        :return: the list of parents of the term
        :rtype: List[GOTerm]
        """
        return [self._ontology.query_term(parent_id) for parent_id in self._parents]

    def get_children(self) -> List["GOTerm"]:
        """get the child terms of the current term

        :return: the list of children of the term
        :rtype: List[GOTerm]
        """
        return [self._ontology.query_term(child_id) for child_id in self._children]


class Neo4jGOOntology(object):
    """ontology with the same properties and methods defined by goatools GODag

    only the properties used for gene descriptions are implemented
    """

    def __init__(self, root_terms: List[Tuple[str, str, List[str]]], node_children_dict, node_parents_dict, node_names,
                 alt_ids):
        self.node_children_dict = node_children_dict
        self.node_parents_dict = node_parents_dict
        self.node_depth = defaultdict(int)
        self.node_names = node_names
        self.alt_ids = alt_ids
        for root_term in root_terms:
            self.calculate_all_depths_in_branch(root_term[0])

    def calculate_all_depths_in_branch(self, root_id: str, current_depth: int = 0):
        """calculate and set depth recursively for all terms in a branch

        :param root_id: the ID of the root term of the branch to process
        :type root_id: str
        :param current_depth: the current depth in the ontology
        :type current_depth: int
        """
        self.node_depth[root_id] = max(self.node_depth[root_id], current_depth)
        for child_id in self.node_children_dict[root_id]:
            self.calculate_all_depths_in_branch(root_id=child_id, current_depth=current_depth + 1)

    def query_term(self, term_id: str):
        """retrieve a term from its ID

        :param term_id: the ID of the term
        :type term_id: str
        :return: the term
        :rtype: GOTerm
        """
        if term_id not in self.node_names:
            if term_id in self.alt_ids:
                term_id = self.alt_ids[term_id]
            else:
                return None
        return GOTerm(name=self.node_names[term_id], depth=self.node_depth[term_id], node_id=term_id,
                      parents=self.node_parents_dict[term_id], children=self.node_children_dict[term_id], ontology=self)


class DataFetcher(metaclass=ABCMeta):
    """retrieve data for gene descriptions from different sources"""

    @abstractmethod
    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str]):
        self.go_data = defaultdict(list)
        self.go_ontology = None
        self.go_terms_exclusion_list = go_terms_exclusion_list
        self.go_terms_replacement_dict = go_terms_replacement_dict

    @abstractmethod
    def get_gene_data(self) -> Gene:
        pass

    @abstractmethod
    def load_go_data(self) -> None:
        pass

    def get_go_annotations(self, geneid: str, include_obsolete: bool = False, include_negative_results: bool = False,
                           priority_list: Iterable = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "IC", "ISS", "ISO",
                                                      "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "IEA"),
                           desc_stats: SingleDescStats = None) -> List[dict]:
        """
        retrieve go annotations for a given gene id and for a given aspect. The annotations are unique for each pair
        <gene_id, go_term_id>. This means that when multiple annotations for the same pair are found in the go data, the
        one with the evidence code with highest priority is returned (see the *priority_list* parameter to set the
        priority according to evidence codes)

        :param geneid: the id of the gene related to the annotations to retrieve, in standard format
        :type geneid: str
        :param include_obsolete: whether to include obsolete annotations
        :type include_obsolete: bool
        :param include_negative_results: whether to include negative results
        :type include_negative_results: bool
        :param priority_list: the priority list for the evidence codes. If multiple annotations with the same go_term
            are found, only the one with highest priority is returned. The first element in the list has the highest
            priority, whereas the last has the lowest. Only annotations with evidence codes in the priority list are
            returned. All other annotations are ignored
        :type priority_list: List[str]
        :param desc_stats: an object containing the description statistics where to save the total number of annotations
            for the gene
        :type desc_stats: SingleDescStats
        :return: the list of go annotations for the given gene
        :rtype: List[GOAnnotation]
        """
        priority_map = dict(zip(priority_list, reversed(range(len(list(priority_list))))))
        annotations = [annotation for annotation in self.go_data[geneid] if (include_obsolete or
                       not annotation["Is_Obsolete"]) and (include_negative_results or "NOT" not in
                                                           annotation["Qualifier"])]
        if desc_stats:
            desc_stats.total_num_go_annotations = len(annotations)
        go_id_selected_annotation = {}
        for annotation in annotations:
            if annotation["Evidence"] in priority_map.keys():
                if annotation["GO_ID"] in go_id_selected_annotation:
                    if priority_map[annotation["Evidence"]] > \
                            priority_map[go_id_selected_annotation[annotation["GO_ID"]]["Evidence"]]:
                        go_id_selected_annotation[annotation["GO_ID"]] = annotation
                else:
                    go_id_selected_annotation[annotation["GO_ID"]] = annotation
        if desc_stats:
            desc_stats.num_prioritized_go_annotations = len(go_id_selected_annotation.keys())
        return [annotation for annotation in go_id_selected_annotation.values()]

    def get_go_ontology(self):
        return self.go_ontology


class RawDataFetcher(DataFetcher):
    """retrieve data for gene descriptions from raw data files"""

    @abstractmethod
    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str],
                 cache_location: str, use_cache: bool = False):
        super().__init__(go_terms_exclusion_list=go_terms_exclusion_list,
                         go_terms_replacement_dict=go_terms_replacement_dict)
        self.chebi_file_url = ""
        self.chebi_file_cache_path = ""
        self.ls_ontology = None
        self.an_ontology = None
        self.gene_data = {}
        self.chebi_ontology = None
        self.use_cache = use_cache
        self.gene_data_cache_path = ""
        self.gene_data_url = ""
        self.go_ontology_cache_path = ""
        self.go_ontology_url = ""
        self.development_ontology_cache_path = ""
        self.development_ontology_url = ""
        self.anatomy_ontology_cache_path = ""
        self.anatomy_ontology_url = ""
        self.go_annotations_cache_path = ""
        self.go_annotations_url = ""
        self.go_id_name = "DB_Object_ID"

    @staticmethod
    def _get_cached_file(cache_path: str, file_source_url):
        if not os.path.isfile(cache_path):
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            urllib.request.urlretrieve(file_source_url, cache_path)
        file_path = cache_path
        if cache_path.endswith(".gz"):
            with gzip.open(cache_path, 'rb') as f_in, open(cache_path.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            file_path = cache_path.replace(".gz", "")
        return file_path

    @abstractmethod
    def _load_gene_data(self) -> None:
        pass

    def get_gene_data(self, include_dead_genes: bool = False, include_pseudo_genes: bool = False) -> Gene:
        """get all gene data from the fetcher, returning one gene per call

        :param include_dead_genes: whether to include dead genes in the results
        :type include_dead_genes: bool
        :param include_pseudo_genes: whether to include pseudo genes in the results
        :type include_dead_genes: bool
        :return: data for one gene per each call, including gene_id and gene_name
        :rtype: Gene
        """
        if len(self.gene_data) == 0:
            self._load_gene_data()
        for gene_id, gene_obj in self.gene_data.items():
            if (include_dead_genes or not gene_obj.dead) and (include_pseudo_genes or not gene_obj.pseudo):
                yield gene_obj

    def _map_ont_term_to_name(self, ont_id):
        t = None
        if ont_id.startswith("WBls:") and self.ls_ontology is not None:
            t = self.ls_ontology.query_term(ont_id)
        elif ont_id.startswith("WBbt:") and self.an_ontology is not None:
            t = self.an_ontology.query_term(ont_id)
        elif ont_id.startswith("GO:") and self.go_ontology is not None:
            t = self.go_ontology.query_term(ont_id)
        elif ont_id.startswith("ChEBI:") and self.chebi_ontology is not None:
            t = self.chebi_ontology.query_term(ont_id)
        elif ont_id.startswith("WB:WBGene"):
            if ont_id in self.gene_data:
                t = self.gene_data[ont_id]
            else:
                return ont_id
        return t.name if t else ont_id

    def load_go_data(self) -> None:
        """read go data and gene ontology. After calling this function, go annotations containing mapped go names can
        be retrieved by using the :meth:`data_fetcher.WBRawDataFetcher.get_go_annotations` function
        """
        goatools.obo_parser.GOTerm.get_parents = get_parents
        self.go_ontology = GODag(self._get_cached_file(file_source_url=self.go_ontology_url,
                                                       cache_path=self.go_ontology_cache_path),
                                 optional_attrs=["relationship"])
        if self.anatomy_ontology_url != "":
            self.an_ontology = GODag(self._get_cached_file(file_source_url=self.anatomy_ontology_url,
                                                           cache_path=self.anatomy_ontology_cache_path),
                                     optional_attrs=["relationship"])
        if self.development_ontology_url != "":
            self.ls_ontology = GODag(self._get_cached_file(file_source_url=self.development_ontology_url,
                                                           cache_path=self.development_ontology_cache_path),
                                     optional_attrs=["relationship"])
        if self.chebi_file_url != "":
            self.chebi_ontology = GODag(self._get_cached_file(file_source_url=self.chebi_file_url,
                                                              cache_path=self.chebi_file_cache_path))
        self._load_gene_data()
        file_path = self._get_cached_file(cache_path=self.go_annotations_cache_path,
                                          file_source_url=self.go_annotations_url)
        lines_to_skip = 0
        with open(file_path) as file:
            while True:
                if file.readline().strip().startswith("!gaf-version:"):
                    break
                lines_to_skip += 1
        with open(file_path) as file:
            for _ in range(lines_to_skip):
                next(file)
            for annotation in gafiterator(file):
                if self.go_ontology.query_term(annotation["GO_ID"]) and \
                        self.go_ontology.query_term(annotation["GO_ID"]).id not in self.go_terms_exclusion_list:
                    mapped_annotation = annotation
                    mapped_annotation["GO_Name"] = self.go_ontology.query_term(mapped_annotation["GO_ID"]).name
                    mapped_annotation["GO_ID"] = self.go_ontology.query_term(mapped_annotation["GO_ID"]).id
                    for regex_to_substitute, regex_target in self.go_terms_replacement_dict.items():
                        mapped_annotation["GO_Name"] = re.sub(regex_to_substitute, regex_target,
                                                              mapped_annotation["GO_Name"])
                    mapped_annotation["Is_Obsolete"] = \
                        self.go_ontology.query_term(mapped_annotation["GO_ID"]).is_obsolete
                    if annotation["Annotation_Extension"] != "":
                        matches = re.findall('(\([^\)]+\))', mapped_annotation["Annotation_Extension"])
                        for match in matches:
                            ext_id = match[1:-1]
                            ext_translation = self._map_ont_term_to_name(ext_id)
                            mapped_annotation["Annotation_Extension"] = mapped_annotation["Annotation_Extension"]\
                                .replace(match, " " + ext_translation)
                        logging.debug(
                            "Found GO annotation with Annotation_Extension: " + str(mapped_annotation))
                    self.go_data[annotation[self.go_id_name]].append(mapped_annotation)


class WBRawDataFetcher(RawDataFetcher):
    """data fetcher for WormBase raw files for a single species"""

    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str],
                 raw_files_source: str, cache_location: str, release_version: str, species: str, project_id: str,
                 use_cache: bool = False, chebi_file_url: str = ""):
        """create a new data fetcher

        :param go_terms_exclusion_list: list of go ids for terms to exclude
        :type go_terms_exclusion_list: List[str]
        :param go_terms_replacement_dict: dictionary to map go terms to be renamed. Term names can be regex
        :type go_terms_replacement_dict: Dict[str, str]
        :param raw_files_source: base url where to fetch the raw files
        :type raw_files_source: str
        :param cache_location: path to cache directory
        :type cache_location: str
        :param release_version: WormBase release version for the input files
        :type release_version: str
        :param species: WormBase species to fetch
        :type species: str
        :param project_id: project id associated with the species
        :type project_id: str
        :param use_cache: whether to use cached files. If cache is empty, files are downloading from source and stored
            in cache
        :type use_cache: bool
        :param chebi_file_url: url where to fetch the chebi file
        :type chebi_file_url: str
        """
        super().__init__(go_terms_exclusion_list=go_terms_exclusion_list,
                         go_terms_replacement_dict=go_terms_replacement_dict, use_cache=use_cache,
                         cache_location=cache_location)
        self.gene_data_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                 project_id, "annotation", species + '.' + project_id +
                                                 '.' + release_version + ".geneIDs.txt.gz")
        self.gene_data_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + project_id + \
                             '/annotation/' + species + '.' + project_id + '.' + release_version + '.geneIDs.txt.gz'
        self.go_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                   "gene_ontology." + release_version + ".obo")
        self.go_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/gene_ontology.' + \
                               release_version + '.obo'
        self.development_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                            "development_ontology." + release_version + ".obo")
        self.development_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/development_ontology.' + \
                                        release_version + '.obo'
        self.anatomy_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                        "anatomy_ontology." + release_version + ".obo")
        self.anatomy_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/anatomy_ontology.' + \
                                    release_version + '.obo'
        self.chebi_file_cache_path = os.path.join(cache_location, "wormbase", "CHEBI", "chebi_lite.obo.gz")
        self.chebi_file_url = chebi_file_url
        self.go_annotations_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                      project_id, "annotation", species + '.' + project_id + '.' +
                                                      release_version + ".go_annotations.gaf.gz")
        self.go_annotations_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + \
                                  project_id + '/annotation/' + species + '.' + project_id + '.' + release_version + \
                                  '.go_annotations.gaf.gz'

    def _load_gene_data(self) -> None:
        """load all gene data"""
        file_path = self._get_cached_file(cache_path=self.gene_data_cache_path, file_source_url=self.gene_data_url)
        with open(file_path) as file:
            for line in file:
                fields = line.strip().split(',')
                name = fields[2] if fields[2] != '' else fields[3]
                self.gene_data[fields[1]] = Gene(fields[1], name, fields[4] == "Dead", False)


class AGRRawDataFetcher(RawDataFetcher):
    """data fetcher for AGR raw files for a single species"""

    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str],
                 raw_files_source: str, cache_location: str, release_version: str, main_file_name: str,
                 bgi_file_name: str, go_annotations_file_name: str, organism_name: str, use_cache: bool = False,
                 chebi_file_url: str = ""):
        """create a new data fetcher

        :param go_terms_exclusion_list: list of go ids for terms to exclude
        :type go_terms_exclusion_list: List[str]
        :param go_terms_replacement_dict: dictionary to map go terms to be renamed. Term names can be regex
        :type go_terms_replacement_dict: Dict[str, str]
        :param raw_files_source: base url where to fetch the raw files
        :type raw_files_source: str
        :param cache_location: path to cache directory
        :type cache_location: str
        :param release_version: WormBase release version for the input files
        :type release_version: str
        :param main_file_name: file name of the main tar.gz file containing gene information
        :type main_file_name: str
        :param bgi_file_name: file name of the bgi file containing gene information
        :type bgi_file_name: str
        :param go_annotations_file_name: file name of the obo file containing go annotations
        :type go_annotations_file_name: str
        :param organism_name: name of the organism
        :type organism_name: str
        :param use_cache: whether to use cached files. If cache is empty, files are downloading from source and stored
            in cache
        :type use_cache: bool
        :param chebi_file_url: url where to fetch the chebi file
        :type chebi_file_url: str
        """
        super().__init__(go_terms_exclusion_list=go_terms_exclusion_list,
                         go_terms_replacement_dict=go_terms_replacement_dict, use_cache=use_cache,
                         cache_location=cache_location)
        self.main_data_cache_path = os.path.join(cache_location, "agr", release_version, "main", main_file_name)
        self.main_data_url = raw_files_source + '/' + main_file_name
        self.bgi_file_name = bgi_file_name
        self.go_ontology_cache_path = os.path.join(cache_location, "agr", release_version, "GO", "go.obo")
        self.go_ontology_url = raw_files_source + '/' + release_version + '/GO/' + 'go.obo'
        self.chebi_file_cache_path = os.path.join(cache_location, "agr", "CHEBI", "chebi_lite.obo.gz")
        self.chebi_file_url = chebi_file_url
        self.go_annotations_cache_path = os.path.join(cache_location, "agr", release_version, "GO", "ANNOT",
                                                      go_annotations_file_name)
        self.go_annotations_url = raw_files_source + '/' + release_version + '/GO/ANNOT/' + go_annotations_file_name
        self.go_id_name = "DB_Object_Symbol"

    def _load_gene_data(self) -> None:
        if not os.path.isfile(self.main_data_cache_path):
            os.makedirs(os.path.dirname(self.main_data_cache_path), exist_ok=True)
            urllib.request.urlretrieve(self.main_data_url, self.main_data_cache_path)
        if not os.path.isfile(os.path.join(os.path.dirname(self.main_data_cache_path), self.bgi_file_name)):
            tar = tarfile.open(self.main_data_cache_path)
            tar.extractall(path=os.path.dirname(self.main_data_cache_path))
        with open(os.path.join(os.path.dirname(self.main_data_cache_path), self.bgi_file_name)) as fileopen:
            bgi_content = json.load(fileopen)
            for gene in bgi_content["data"]:
                self.gene_data[gene["symbol"]] = Gene(gene["symbol"], gene["symbol"], False, False)


class AGRDBDataFetcher(DataFetcher):
    """data fetcher for AGR neo4j database for a single species"""

    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str],
                 data_provider: str, db_graph, go_ontology=None):
        super().__init__(go_terms_exclusion_list=go_terms_exclusion_list,
                         go_terms_replacement_dict=go_terms_replacement_dict)
        self.db_graph = db_graph
        self.data_provider = data_provider
        self.go_ontology = go_ontology

    @staticmethod
    def query_db(db_graph, query: str, parameters: Dict = None):
        """query the neo4j db

        :param db_graph: a neo4j graph database
        :param query: a cypher
        :type query: str
        :param parameters: a dictionary of the parameters of the query
        :type parameters: Dict
        """
        with db_graph.session() as session:
            with session.begin_transaction() as tx:
                return_set = tx.run(query, parameters=parameters)
        return return_set

    def get_gene_data(self):
        """get all gene data from the fetcher, returning one gene per call

        while reading data for a gene, all related annotations are loaded from neo4j into memory

        :return: data for one gene per each call, including gene_id and gene_name
        :rtype: Gene
        """
        db_query = "match (g:Gene) where g.dataProvider = {dataProvider} return g.symbol, g.primaryKey"
        result_set = self.query_db(db_graph=self.db_graph, query=db_query,
                                   parameters={"dataProvider": self.data_provider})
        for result in result_set:
            annot_db_query = "match (g:Gene)-[annot:ANNOTATED_TO]->(go_term:GOTerm:Ontology) " \
                             "where g.primaryKey = {gene_id} " \
                             "return g.primaryKey, g.symbol, annot.evidence_code, annot.aspect, annot.qualifier, " \
                             "go_term.primaryKey, go_term.name, go_term.is_obsolete"
            annot_result_set = self.query_db(db_graph=self.db_graph, query=annot_db_query,
                                             parameters={"gene_id": result["g.primaryKey"]})
            for annot in annot_result_set:
                if not annot["go_term.is_obsolete"] == "true":
                    onto_node = self.go_ontology.query_term(annot["go_term.primaryKey"])
                    if onto_node:
                        self.go_data[result["g.symbol"]].append({
                            "DB": self.data_provider,
                            "DB_Object_ID": annot["g.primaryKey"],
                            "DB_Object_Symbol": annot["g.symbol"],
                            "Qualifier": annot["annot.qualifier"],
                            "GO_ID": onto_node.id,
                            "DB:Reference": None,
                            "Evidence": annot["annot.evidence_code"],
                            "With": "",
                            "Aspect": annot["annot.aspect"],
                            "DB_Object_Name": None,
                            "Synonym": None,
                            "DB_Object_Type": None,
                            "Taxon_ID": None,
                            "Date": None,
                            "Assigned_By": None,
                            "Annotation_Extension": None,
                            "Gene_Product_Form_ID": None,
                            "GO_Name": onto_node.name,
                            "Is_Obsolete": False
                        })
                    else:
                        logging.warning("Annotated GO Term " + annot["go_term.primaryKey"] + " not found in GO ontology")
            yield Gene(result["g.primaryKey"], result["g.symbol"], False, False)

    def load_go_data(self):
        """load GO ontology from neo4j if not already loaded"""
        if not self.go_ontology:

            # get root terms
            db_query = "match path=(child:GOTerm:Ontology)-[r:IS_A|PART_OF]->(parent:GOTerm:Ontology) " \
                       "where size((parent)-[:IS_A|PART_OF]->()) = 0 " \
                       "return distinct parent.primaryKey, parent.alt_ids, parent.name"
            result_set = self.query_db(db_graph=self.db_graph, query=db_query)
            root_terms = []
            for result in result_set:
                root_terms.append((result["parent.primaryKey"], result["parent.name"], result["parent.alt_ids"]))

            # cypher query to build go_ontology object
            db_query = "match (child:GOTerm:Ontology)-[:IS_A|PART_OF]->(parent:GOTerm:Ontology) " \
                       "return distinct child.primaryKey, child.alt_ids, child.name, parent.primaryKey"
            result_set = self.query_db(db_graph=self.db_graph, query=db_query)
            node_parents_dict = defaultdict(list)
            node_children_dict = defaultdict(list)
            node_names = {}
            alt_ids = {}
            for result in result_set:
                node_parents_dict[result["child.primaryKey"]].append(result["parent.primaryKey"])
                node_children_dict[result["parent.primaryKey"]].append(result["child.primaryKey"])
                node_names[result["child.primaryKey"]] = result["child.name"]
                for alt_id in list(result["child.alt_ids"]):
                    alt_ids[alt_id] = result["child.primaryKey"]
            for root_term_id, root_term_name, root_term_alt_ids in root_terms:
                node_names[root_term_id] = root_term_name
                for alt_id in list(root_term_alt_ids):
                    alt_ids[alt_id] = root_term_id

            self.go_ontology = Neo4jGOOntology(root_terms=root_terms, node_children_dict=node_children_dict,
                                               node_parents_dict=node_parents_dict, node_names=node_names,
                                               alt_ids=alt_ids)
