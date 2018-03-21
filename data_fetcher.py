import gzip
import json
import tarfile
import urllib.request
import shutil
import os
import logging
import re

from goatools.obo_parser import GODag
from Bio.UniProt.GOA import gafiterator
from abc import ABCMeta, abstractmethod
from collections import namedtuple, defaultdict
from typing import List, Iterable, Dict

from descriptions_writer import SingleDescStats

Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


class RawDataFetcher(metaclass=ABCMeta):

    @abstractmethod
    def __init__(self, go_terms_exclusion_list: List[str], go_terms_replacement_dict: Dict[str, str],
                 cache_location: str, use_cache: bool = False):
        self.chebi_file_url = ""
        self.chebi_file_cache_path = ""
        self.go_data = defaultdict(list)
        self.go_ontology = None
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
        self.go_terms_exclusion_list = go_terms_exclusion_list
        self.go_terms_replacement_dict = go_terms_replacement_dict
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
        self.go_ontology = GODag(self._get_cached_file(file_source_url=self.go_ontology_url,
                                                       cache_path=self.go_ontology_cache_path))
        if self.anatomy_ontology_url != "":
            self.an_ontology = GODag(self._get_cached_file(file_source_url=self.anatomy_ontology_url,
                                                           cache_path=self.anatomy_ontology_cache_path))
        if self.development_ontology_url != "":
            self.ls_ontology = GODag(self._get_cached_file(file_source_url=self.development_ontology_url,
                                                           cache_path=self.development_ontology_cache_path))
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



