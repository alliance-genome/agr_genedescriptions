import gzip
import urllib.request
from collections import namedtuple, defaultdict
from enum import Enum
import os
from urllib.parse import urlparse
from typing import List
import logging
import re

Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])
GOAnnotation = namedtuple('GOAnnotation', ['go_id', 'qualifier', 'paper_reference', 'evidence_code', 'aspect',
                                           'annotation_ext', 'go_name', 'is_obsolete'])


class GO_ASPECT(Enum):
    MOLECULAR_FUNCTION = 0
    BIOLOGICAL_PROCESS = 1
    CELLULAR_COMPONENT = 2


class WBRawDataFetcher:
    """data fetcher for WormBase raw files for a single species"""
    def __init__(self, raw_files_source: str, chebi_files_source: str, cache_location: str, release_version: str,
                 species: str, project_id: str, use_cache: bool, ):
        """create a new data fetcher

        :param raw_files_source: base url where to fetch the raw files
        :type raw_files_source: str
        :param chebi_files_source: base url where to fetch the chebi files
        :type chebi_files_source: str
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
        """
        self.raw_files_source = raw_files_source
        self.chebi_files_source = chebi_files_source
        self.cache_location = cache_location
        self.release_version = release_version
        self.species = species
        self.project_id = project_id
        self.go_data = defaultdict(set)
        self.go_ontology = {}
        self.ls_ontology = {}
        self.an_ontology = {}
        self.gene_data = {}
        self.chebi_ontology = {}
        self.use_cache = use_cache

    def _fill_cache_if_empty_and_activated(self, cache_url, file_source_url):
        cache_url_parsed = urlparse(cache_url)
        if self.use_cache and not os.path.isfile(cache_url_parsed.path):
            os.makedirs(os.path.dirname(cache_url_parsed.path), exist_ok=True)
            urllib.request.urlretrieve(file_source_url, cache_url_parsed.path)

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

    def _load_gene_data(self) -> None:
        """load all gene data"""
        cache_url = os.path.join(self.cache_location, self.release_version, "species", self.species, self.project_id,
                                 "annotation", self.species + '.' + self.project_id + '.' + self.release_version +
                                 ".geneIDs.txt.gz")
        source_address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                         self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + \
                         self.release_version + '.geneIDs.txt.gz'
        self._fill_cache_if_empty_and_activated(cache_url=cache_url, file_source_url=source_address)
        address = cache_url if self.use_cache else source_address
        with urllib.request.urlopen(address) as url:
            gzip_file = gzip.GzipFile(fileobj=url)
            for line in gzip_file:
                fields = line.decode("utf-8").strip().split(',')
                name = fields[2] if fields[2] != '' else fields[3]
                self.gene_data[fields[1]] = Gene(fields[1], name, fields[4] == "Dead", False)

    def load_go_data(self) -> None:
        """read go data and gene ontology. After calling this function, go annotations containing mapped go names can
        be retrieved by using the :meth:`data_fetcher.WBRawDataFetcher.get_go_annotations` function
        """
        self._load_ontology_data(url=self.raw_files_source + '/' + self.release_version + '/ONTOLOGY/gene_ontology.' + \
                                 self.release_version + '.obo', cache_path=os.path.join(self.cache_location,
                                                                                        self.release_version,
                                                                                        "ONTOLOGY", "gene_ontology." +
                                                                                        self.release_version + ".obo"),
                                 field_names=["name", "is_obsolete"], ontology_dict=self.go_ontology, gzip_file=False)
        self._load_ontology_data(url=self.raw_files_source + '/' + self.release_version +
                                     '/ONTOLOGY/development_ontology.' + self.release_version + '.obo',
                                 cache_path=os.path.join(self.cache_location, self.release_version, "ONTOLOGY",
                                                         "development_ontology." + self.release_version + ".obo"),
                                 field_names=["name"], ontology_dict=self.ls_ontology, gzip_file=False)
        self._load_ontology_data(url=self.raw_files_source + '/' + self.release_version +
                                     '/ONTOLOGY/anatomy_ontology.' + self.release_version + '.obo',
                                 cache_path=os.path.join(self.cache_location, self.release_version, "ONTOLOGY",
                                                         "anatomy_ontology." + self.release_version + ".obo"),
                                 field_names=["name"], ontology_dict=self.an_ontology, gzip_file=False)
        self._load_ontology_data(url=self.chebi_files_source + '/' + 'chebi_lite.obo.gz',
                                 cache_path=os.path.join(self.cache_location, "CHEBI", "chebi_lite.obo.gz"),
                                 field_names=["name"], ontology_dict=self.chebi_ontology, gzip_file=True)
        self._load_gene_data()
        cache_url = os.path.join(self.cache_location, self.release_version, "species", self.species, self.project_id,
                                 "annotation", self.species + '.' + self.project_id + '.' + self.release_version +
                                 ".go_annotations.gaf.gz")
        source_address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                         self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + \
                         self.release_version + '.go_annotations.gaf.gz'
        self._fill_cache_if_empty_and_activated(cache_url=cache_url, file_source_url=source_address)
        address = cache_url if self.use_cache else source_address
        with urllib.request.urlopen(address) as url:
            gzip_file = gzip.GzipFile(fileobj=url)
            for line in gzip_file:
                line = line.decode("utf-8")
                if not line.startswith("!"):
                    fields = line.strip("\n").split('\t')
                    go_aspect = None
                    if fields[8] == 'C':
                        go_aspect = GO_ASPECT.CELLULAR_COMPONENT
                    elif fields[8] == 'F':
                        go_aspect = GO_ASPECT.MOLECULAR_FUNCTION
                    elif fields[8] == 'P':
                        go_aspect = GO_ASPECT.BIOLOGICAL_PROCESS
                    is_obsolete = False
                    if "is_obsolete" in self.go_ontology[fields[4]]["is_obsolete"] and \
                            self.go_ontology[fields[4]]["is_obsolete"] == 'true':
                        is_obsolete = True
                    extension = ""
                    if fields[15] != "":
                        extension = fields[15]
                        matches = re.findall('(\([^\)]+\))', fields[15])
                        for match in matches:
                            ext_id = match[1:-1]
                            ext_translation = self._ontology_translation(ext_id)
                            extension = extension.replace(match, " " + ext_translation)
                        logging.debug("Found GO annotation with field 16: " + fields[1] + "\t" + fields[2] + "\t" +
                                      fields[3] + "\t" + self._map_id_ontology_name(fields[4], self.go_ontology) +
                                      "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8] + "\t" + fields[14] +
                                      "\t" + extension)
                    self.go_data[fields[1]].add(GOAnnotation(fields[4], fields[3], fields[5], fields[6], go_aspect,
                                                             extension, self.go_ontology[fields[4]]["name"],
                                                             is_obsolete))

    @staticmethod
    def _map_id_ontology_name(ont_id: str, ontology_dict: dict):
        if ont_id in ontology_dict:
            return ontology_dict[ont_id]["name"]
        else:
            return ont_id

    def _ontology_translation(self, ont_id):
        if ont_id.startswith("WBls:"):
            return self._map_id_ontology_name(ont_id, self.ls_ontology)
        elif ont_id.startswith("WBbt:"):
            return self._map_id_ontology_name(ont_id, self.an_ontology)
        elif ont_id.startswith("GO:"):
            return self._map_id_ontology_name(ont_id, self.go_ontology)
        elif ont_id.startswith("ChEBI:"):
            return self._map_id_ontology_name(ont_id, self.chebi_ontology)
        elif ont_id.startswith("WB:WBGene"):
            if ont_id in self.gene_data:
                return self.gene_data[ont_id].name
            else:
                return ont_id
        else:
            return ont_id

    def _load_ontology_data(self, url: str, cache_path: str, field_names: List[str], ontology_dict: dict,
                            gzip_file: bool = False):
        """read ontology data"""
        self._fill_cache_if_empty_and_activated(cache_url=cache_path, file_source_url=url)
        address = cache_path if self.use_cache else url
        with urllib.request.urlopen(address) as url:
            if gzip_file:
                url = gzip.GzipFile(fileobj=url)
            ont_id = None
            ont_entry = defaultdict(str)
            for line in url:
                line = line.decode("utf-8")
                if line.strip() == '[Term]':
                    if ont_id is not None:
                        ontology_dict[ont_id] = ont_entry
                        ont_id = None
                        ont_entry = {"is_obsolete": "false"}
                else:
                    fields = line.strip().split(": ")
                    if fields[0] == "id":
                        ont_id = fields[1]
                    else:
                        for fn in field_names:
                            if fields[0] == fn:
                                ont_entry[fn] = fields[1]

    def get_go_annotations(self, geneid: str, include_obsolete: bool = False,
                           priority_list: tuple = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "IC", "ISS", "ISO", "ISA",
                                                   "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "IEA"
                                                   )) -> List[GOAnnotation]:
        """retrieve go annotations for a given gene id and for a given aspect. The annotations are unique for each pair
        <gene_id, go_term_id>. This means that when multiple annotations for the same pair are found in the go data, the
        one with the evidence code with highest priority is returned (see the *priority_list* parameter to set the
        priority according to evidence codes)

        :param geneid: the id of the gene related to the annotations to retrieve, in standard format
        :type geneid: str
        :param include_obsolete: whether to include obsolete annotations
        :type include_obsolete: bool
        :param priority_list: the priority list for the evidence codes. If multiple annotations with the same go_term
            are found, only the one with highest priority is returned. The first element in the list has the highest
            priority, whereas the last has the lowest. Only annotations with evidence codes in the priority list are
            returned. All other annotations are ignored
        :type priority_list: List[str]
        :return: the list of go annotations for the given gene
        :rtype: List[GOAnnotation]
        """
        priority_map = dict(zip(priority_list, reversed(range(len(priority_list)))))
        annotations = [annotation for annotation in self.go_data[geneid] if include_obsolete or
                       not annotation.is_obsolete]
        go_id_selected_annotation = {}
        for annotation in annotations:
            if annotation.evidence_code in priority_map.keys():
                if annotation.go_id in go_id_selected_annotation:
                    if priority_map[annotation.evidence_code] > \
                            priority_map[go_id_selected_annotation[annotation.go_id].evidence_code]:
                        go_id_selected_annotation[annotation.go_id] = annotation
                else:
                    go_id_selected_annotation[annotation.go_id] = annotation

        return [annotation for annotation in go_id_selected_annotation.values()]

