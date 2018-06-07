import gzip
import json
import tarfile
import urllib.request
import shutil
import os
import re
from enum import Enum
from abc import ABCMeta, abstractmethod
from collections import namedtuple
from typing import List, Iterable, Dict
from ontobio import AssociationSetFactory
from genedescriptions.descriptions_rules import SingleDescStats, set_all_depths_in_subgraph
from ontobio.ontol_factory import OntologyFactory
from ontobio.ontol import Ontology
from ontobio.assocmodel import AssociationSet
from ontobio.io.gafparser import GafParser

Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


class DataType(Enum):
    GO = 1
    DO = 2


class DataFetcher(metaclass=ABCMeta):
    """retrieve data for gene descriptions from different sources"""

    @abstractmethod
    def __init__(self, go_relations: List[str] = None, do_relations: List[str] = None, use_cache: bool = False):
        self.go_associations: AssociationSet = None
        self.go_ontology: Ontology = None
        self.do_ontology: Ontology = None
        self.do_associations: AssociationSet = None
        self.gene_data: Dict[str, Gene] = None
        self.go_relations = go_relations
        self.do_relations = do_relations
        self.use_cache = use_cache

    def _get_cached_file(self, cache_path: str, file_source_url):
        if not os.path.isfile(cache_path):
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            urllib.request.urlretrieve(file_source_url, cache_path)
        elif not self.use_cache:
            urllib.request.urlretrieve(file_source_url, cache_path)
        file_path = cache_path
        if cache_path.endswith(".gz"):
            with gzip.open(cache_path, 'rb') as f_in, open(cache_path.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            file_path = cache_path.replace(".gz", "")
        return file_path

    def get_gene_data(self, include_dead_genes: bool = False, include_pseudo_genes: bool = False) -> Gene:
        """get all gene data from the fetcher, returning one gene per call

        :param include_dead_genes: whether to include dead genes in the results
        :type include_dead_genes: bool
        :param include_pseudo_genes: whether to include pseudo genes in the results
        :type include_dead_genes: bool
        :return: data for one gene per each call, including gene_id and gene_name
        :rtype: Gene
        """
        if self.gene_data and len(self.gene_data) > 0:
            for gene_id, gene_obj in self.gene_data.items():
                if (include_dead_genes or not gene_obj.dead) and (include_pseudo_genes or not gene_obj.pseudo):
                    yield gene_obj

    @staticmethod
    def remove_blacklisted_annotations(association_set: AssociationSet, ontology: Ontology,
                                       terms_blacklist) -> AssociationSet:
        """remove annotations linked to blacklisted ontology terms from an association set

        :param association_set: the original association set
        :type association_set: AssociationSet
        :param ontology: the ontology linked to the annotations
        :type ontology: Ontology
        :param terms_blacklist: the list of ontology terms related to the annotations to be removed
        :type terms_blacklist: List[str]
        :return: the filtered annotations
        :rtype: AssociationSet
        """
        associations = []
        for subj_associations in association_set.associations_by_subj.values():
            for association in subj_associations:
                if association["object"]["id"] not in terms_blacklist:
                    associations.append(association)
        return AssociationSetFactory().create_from_assocs(assocs=associations, ontology=ontology)

    @staticmethod
    def rename_ontology_terms(ontology: Ontology, terms_replacement_regex: Dict[str, str]) -> None:
        """rename ontology terms based on regular expression matching

        :param ontology: the ontology containing the terms to be renamed
        :type ontology: Ontology
        :param terms_replacement_regex: a dictionary containing the regular expression to be applied for renaming terms.
            Each key must be a regular expression to search for terms and the associated value another regular
            expression that defines the final result
        :type terms_replacement_regex: Dict[str, str]
        """
        for regex_to_substitute, regex_target in terms_replacement_regex.items():
            for node in ontology.search(regex_to_substitute, is_regex=True):
                ontology.node(node)["label"] = re.sub(regex_to_substitute, regex_target, ontology.node(node)["label"])

    def set_ontology(self, ontology_type: DataType, ontology: Ontology,
                     terms_replacement_regex: Dict[str, str]) -> None:
        """set the go ontology and apply terms renaming

        :param ontology_type: the type of ontology to set
        :type ontology_type: DataType
        :param ontology: an ontology object to set as go ontology
        :type ontology: Ontology
        :param terms_replacement_regex: a dictionary containing the regular expression to be applied for renaming terms.
            Each key must be a regular expression to search for terms and the associated value another regular
            expression that defines the final result
        :type terms_replacement_regex: Dict[str, str]
        """
        new_ontology = None
        if ontology_type == DataType.GO:
            self.go_ontology = ontology.subontology(relations=self.go_relations)
            new_ontology = self.go_ontology
        elif ontology_type == DataType.DO:
            self.do_ontology = ontology.subontology(relations=self.do_relations)
            new_ontology = self.do_ontology
        self.rename_ontology_terms(ontology=new_ontology, terms_replacement_regex=terms_replacement_regex)
        for root_id in new_ontology.get_roots():
            set_all_depths_in_subgraph(ontology=new_ontology, root_id=root_id, relations=None)

    def load_ontology_from_file(self, ontology_type: DataType, ontology_url: str, ontology_cache_path: str,
                                terms_replacement_regex: Dict[str, str] = None) -> None:
        """load go ontology from file

        :param ontology_type: the type of ontology to set
        :type ontology_type: DataType
        :param ontology_url: url to the ontology file
        :type ontology_url: str
        :param ontology_cache_path: path to cache file for the ontology
        :type ontology_cache_path: str
        :param terms_replacement_regex: a dictionary containing the regular expression to be applied for renaming terms.
            Each key must be a regular expression to search for terms and the associated value another regular
            expression that defines the final result
        :type terms_replacement_regex: Dict[str, str]
        """
        new_ontology = None
        if ontology_type == DataType.GO:
            self.go_ontology = OntologyFactory().create(self._get_cached_file(file_source_url=ontology_url,
                                                                              cache_path=ontology_cache_path)
                                                        ).subontology(relations=self.go_relations)
            new_ontology = self.go_ontology
        elif ontology_type == DataType.DO:
            self.do_ontology = OntologyFactory().create(self._get_cached_file(file_source_url=ontology_url,
                                                                              cache_path=ontology_cache_path)
                                                        ).subontology(relations=self.do_relations)
            new_ontology = self.do_ontology
        if terms_replacement_regex:
            self.rename_ontology_terms(ontology=new_ontology, terms_replacement_regex=terms_replacement_regex)
        for root_id in new_ontology.get_roots():
            set_all_depths_in_subgraph(ontology=new_ontology, root_id=root_id, relations=None)

    def set_associations(self, associations_type: DataType, associations: AssociationSet,
                         exclusion_list: List[str]) -> None:
        """set the go annotations and remove blacklisted annotations

        :param associations_type: the type of associations to set
        :type associations_type: DataType
        :param associations: an association object to set as go annotations
        :type associations: AssociationSet
        :param exclusion_list: the list of ontology terms related to the annotations to be removed
        :type exclusion_list: List[str]
        """
        if associations_type == DataType.GO:
            self.go_associations = self.remove_blacklisted_annotations(association_set=associations,
                                                                       ontology=self.go_ontology,
                                                                       terms_blacklist=exclusion_list)
        elif associations_type == DataType.DO:
            self.do_associations = self.remove_blacklisted_annotations(association_set=associations,
                                                                       ontology=self.do_ontology,
                                                                       terms_blacklist=exclusion_list)

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, exclusion_list: List[str]) -> None:
        """load go associations from file

        :param associations_type: the type of associations to set
        :type associations_type: DataType
        :param associations_url: url to the association file
        :type associations_url: str
        :param associations_cache_path: path to cache file for the associations
        :type associations_cache_path: str
        :param exclusion_list: the list of ontology terms related to the annotations to be removed
        :type exclusion_list: List[str]
        """
        if associations_type == DataType.GO:
            self.go_associations = AssociationSetFactory().create_from_assocs(assocs=GafParser().parse(
                file=self._get_cached_file(cache_path=associations_cache_path, file_source_url=associations_url),
                skipheader=True), ontology=self.go_ontology)
            self.go_associations = self.remove_blacklisted_annotations(association_set=self.go_associations,
                                                                       ontology=self.go_ontology,
                                                                       terms_blacklist=exclusion_list)
        elif associations_type == DataType.DO:
            self.do_associations = AssociationSetFactory().create_from_assocs(assocs=GafParser().parse(
                file=self._get_cached_file(cache_path=associations_cache_path, file_source_url=associations_url),
                skipheader=True), ontology=self.do_ontology)
            self.do_associations = self.remove_blacklisted_annotations(association_set=self.do_associations,
                                                                       ontology=self.do_ontology,
                                                                       terms_blacklist=exclusion_list)

    def get_annotations_for_gene(self, gene_id: str, annot_type: DataType = DataType.GO,
                                 include_obsolete: bool = False, include_negative_results: bool = False,
                                 priority_list: Iterable = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "IC", "ISS",
                                                            "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD",
                                                            "RCA", "IEA"),
                                 desc_stats: SingleDescStats = None) -> List[Dict]:
        """
        retrieve go annotations for a given gene id and a given type. The annotations are unique for each pair
        <gene_id, go_term_id>. This means that when multiple annotations for the same pair are found in the go data, the
        one with the evidence code with highest priority is returned (see the *priority_list* parameter to set the
        priority according to evidence codes)

        :param gene_id: the id of the gene related to the annotations to retrieve, in standard format
        :type gene_id: str
        :param annot_type: type of annotations to read
        :type annot_type: DataType
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
        :rtype: List[Dict]
        """
        dataset = None
        ontology = None
        if annot_type == DataType.GO:
            dataset = self.go_associations
            ontology = self.go_ontology
        elif annot_type == DataType.DO:
            dataset = self.do_associations
            ontology = self.do_ontology
        priority_map = dict(zip(priority_list, reversed(range(len(list(priority_list))))))
        annotations = [annotation for annotation in dataset.associations(gene_id) if (include_obsolete or
                                                                                      not ontology.is_obsolete(
                                                                                          annotation["object"]["id"]))
                       and (include_negative_results or "NOT" not in annotation["qualifiers"])]
        if desc_stats:
            desc_stats.total_num_go_annotations = len(annotations)
        id_selected_annotation = {}
        for annotation in annotations:
            if annotation["evidence"]["type"] in priority_map.keys():
                if annotation["object"]["id"] in id_selected_annotation:
                    if priority_map[annotation["evidence"]["type"]] > \
                            priority_map[id_selected_annotation[annotation["object"]["id"]]["evidence"]["type"]]:
                        id_selected_annotation[annotation["object"]["id"]] = annotation
                else:
                    id_selected_annotation[annotation["object"]["id"]] = annotation
        if desc_stats:
            desc_stats.num_prioritized_go_annotations = len(id_selected_annotation.keys())
        return [annotation for annotation in id_selected_annotation.values()]

    def set_gene_data(self, gene_data: List[Gene]):
        for gene in gene_data:
            self.gene_data[gene.id] = gene

    @abstractmethod
    def load_gene_data_from_file(self):
        pass


class WBDataFetcher(DataFetcher):
    """data fetcher for WormBase raw files for a single species"""

    def __init__(self, raw_files_source: str, cache_location: str, release_version: str, species: str, project_id: str,
                 go_relations: List[str] = None, do_relations: List[str] = None, use_cache: bool = False):
        """create a new data fetcher for WormBase. Files will be downloaded from WB ftp site. For convenience, file
        locations are automatically generated and stored in class variables ending in _url for remote filed and
        _cache_path for caching

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
        """
        super().__init__(go_relations=go_relations, do_relations=do_relations, use_cache=use_cache)
        self.gene_data_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                 project_id, "annotation", species + '.' + project_id +
                                                 '.' + release_version + ".geneIDs.txt.gz")
        self.gene_data_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + project_id + \
                             '/annotation/' + species + '.' + project_id + '.' + release_version + '.geneIDs.txt.gz'
        self.go_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                   "gene_ontology." + release_version + ".obo")
        self.go_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/gene_ontology.' + \
                               release_version + '.obo'
        self.go_associations_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                       project_id, "annotation", species + '.' + project_id + '.' +
                                                       release_version + ".go_annotations.gaf.gz")
        self.go_associations_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + \
                                   project_id + '/annotation/' + species + '.' + project_id + '.' + release_version + \
                                  '.go_annotations.gaf.gz'
        self.do_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/disease_ontology.' + \
                               release_version + '.obo'
        self.do_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                   "disease_ontology." + release_version + ".obo")
        self.do_associations_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                       project_id, "annotation", species + '.' + project_id + '.' +
                                                       release_version + ".do_annotations.wb")
        self.do_associations_url = raw_files_source + '/' + release_version + '/ONTOLOGY/disease_association.' + \
                                   release_version + '.wb'
        self.do_associations_new_cache_path = os.path.join(cache_location, "wormbase", release_version, "species",
                                                           species, project_id, "annotation", species + '.' +
                                                           project_id + '.' + release_version +
                                                          ".do_annotations.daf.txt")
        self.do_associations_new_url = raw_files_source + '/' + release_version + '/ONTOLOGY/disease_association.' + \
                                       release_version + '.daf.txt'

    def load_gene_data_from_file(self) -> None:
        if not self.gene_data or len(self.gene_data.items()) == 0:
            self.gene_data = {}
            file_path = self._get_cached_file(cache_path=self.gene_data_cache_path, file_source_url=self.gene_data_url)
            with open(file_path) as file:
                for line in file:
                    fields = line.strip().split(',')
                    name = fields[2] if fields[2] != '' else fields[3]
                    self.gene_data["WB:" + fields[1]] = Gene("WB:" + fields[1], name, fields[4] == "Dead", False)

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, exclusion_list: List[str]) -> None:
        if associations_type == DataType.GO:
            super().load_associations_from_file(associations_type=associations_type, associations_url=associations_url,
                                                associations_cache_path=associations_cache_path,
                                                exclusion_list=exclusion_list)
        elif associations_type == DataType.DO:
            self.do_associations = AssociationSetFactory().create_from_assocs(
                assocs=GafParser().parse(file=self._get_cached_file(cache_path=self.do_associations_cache_path,
                                                                    file_source_url=self.do_associations_url),
                                         skipheader=True),
                ontology=self.do_ontology)
            associations = []
            for subj_associations in self.do_associations.associations_by_subj.values():
                for association in subj_associations:
                    if association["evidence"]["type"] == "IEA":
                        associations.append(association)
            file_path = self._get_cached_file(cache_path=self.do_associations_new_cache_path,
                                              file_source_url=self.do_associations_new_url)
            header = True
            for line in open(file_path):
                if not line.strip().startswith("!"):
                    if not header:
                        linearr = line.strip().split("\t")
                        if self.do_ontology.node(linearr[10]) and linearr[16] != "IEA":
                            associations.append({"source_line": line,
                                                 "subject": {
                                                     "id": linearr[2],
                                                     "label": linearr[3],
                                                     "type": line[1],
                                                     "fullname": "",
                                                     "synonyms": [],
                                                     "taxon": {"id": linearr[0]}

                                                 },
                                                 "object": {
                                                     "id": linearr[10],
                                                     "taxon": ""
                                                 },
                                                 "qualifiers": linearr[9].split("|"),
                                                 "aspect": "D",
                                                 "relation": {"id": None},
                                                 "negated": False,
                                                 "evidence": {
                                                     "type": linearr[16],
                                                     "has_supporting_reference": linearr[18].split("|"),
                                                     "with_support_from": [],
                                                     "provided_by": linearr[20],
                                                     "date": linearr[19]
                                                     }
                                                 })
                    else:
                        header = False
            self.do_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                              ontology=self.do_ontology)
            self.do_associations = self.remove_blacklisted_annotations(association_set=self.do_associations,
                                                                       ontology=self.do_ontology,
                                                                       terms_blacklist=exclusion_list)

    def load_all_data_from_file(self, go_terms_replacement_regex: Dict[str, str] = None,
                                go_terms_exclusion_list: List[str] = None,
                                do_terms_replacement_regex: Dict[str, str] = None,
                                do_terms_exclusion_list: List[str] = None) -> None:
        """load all data types from stored file locations

        :param go_terms_replacement_regex: a dictionary containing the regular expression to be applied for renaming
            GO terms. Each key must be a regular expression to search for terms and the associated value another regular
            expression that defines the final result
        :type go_terms_replacement_regex: Dict[str, str]
        :param go_terms_exclusion_list: the list of GO ontology terms related to the GO annotations to be removed
        :type go_terms_exclusion_list: List[str]
        :param do_terms_replacement_regex: a dictionary containing the regular expression to be applied for renaming
            DO terms. Each key must be a regular expression to search for terms and the associated value another regular
            expression that defines the final result
        :type do_terms_replacement_regex: Dict[str, str]
        :param do_terms_exclusion_list: the list of DO ontology terms related to the DO annotations to be removed
        :type do_terms_exclusion_list: List[str]
        """
        self.load_gene_data_from_file()
        self.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=self.go_ontology_url,
                                     ontology_cache_path=self.go_ontology_cache_path,
                                     terms_replacement_regex=go_terms_replacement_regex)
        self.load_associations_from_file(associations_type=DataType.GO, associations_url=self.go_associations_url,
                                         associations_cache_path=self.go_associations_cache_path,
                                         exclusion_list=go_terms_exclusion_list)
        self.load_ontology_from_file(ontology_type=DataType.DO, ontology_url=self.do_ontology_url,
                                     ontology_cache_path=self.do_ontology_cache_path,
                                     terms_replacement_regex=do_terms_replacement_regex)
        self.load_associations_from_file(associations_type=DataType.DO, associations_url=self.do_associations_url,
                                         associations_cache_path=self.do_associations_cache_path,
                                         exclusion_list=do_terms_exclusion_list)


class AGRDataFetcher(DataFetcher):
    """data fetcher for AGR raw files for a single species"""

    def __init__(self, raw_files_source: str, cache_location: str, release_version: str, main_file_name: str,
                 bgi_file_name: str, go_annotations_file_name: str, organism_name: str, go_relations: List[str] = None,
                 do_relations: List[str] = None, use_cache: bool = False):
        """create a new data fetcher

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

        """
        super().__init__(do_relations=do_relations, go_relations=go_relations, use_cache=use_cache)
        self.main_data_cache_path = os.path.join(cache_location, "agr", release_version, "main", main_file_name)
        self.main_data_url = raw_files_source + '/' + main_file_name
        self.bgi_file_name = bgi_file_name
        self.go_ontology_cache_path = os.path.join(cache_location, "agr", release_version, "GO", "go.obo")
        self.go_ontology_url = raw_files_source + '/' + release_version + '/GO/' + 'go.obo'
        self.go_annotations_cache_path = os.path.join(cache_location, "agr", release_version, "GO", "ANNOT",
                                                      go_annotations_file_name)
        self.go_annotations_url = raw_files_source + '/' + release_version + '/GO/ANNOT/' + go_annotations_file_name
        #self.go_id_name = "DB_Object_Symbol"

    def load_gene_data_from_file(self) -> None:
        if len(self.gene_data.items()) == 0:
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
