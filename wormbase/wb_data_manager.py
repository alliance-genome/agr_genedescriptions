import os
import inflect

from collections import defaultdict
from typing import List
from ontobio import AssociationSetFactory
from ontobio.io.gafparser import GafParser
from genedescriptions.commons import DataType, Gene, Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import ExpressionClusterFeature, DataManager, ExpressionClusterType


class WBDataManager(DataManager):
    """data fetcher for WormBase raw files for a single species"""

    def __init__(self, raw_files_source: str, cache_location: str, release_version: str, species: str, project_id: str,
                 go_relations: List[str] = None, do_relations: List[str] = None, use_cache: bool = False,
                 sister_sp_fullname: str = "", human_orthologs_go_ontology_url: str = None,
                 human_orthologs_go_associations_url: str = None, expression_cluster_anatomy_prefix: str = None,
                 expression_cluster_molreg_prefix: str = None):
        """create a new data fetcher for WormBase. Files will be downloaded from WB ftp site. For convenience, file
        locations are automatically generated and stored in class variables ending in _url for remote filed and
        _cache_path for caching

        Args:
            raw_files_source (str): base url where to fetch the raw files
            cache_location (str): path to cache directory
            release_version (str): WormBase release version for the input files
            species (str): WormBase species to fetch
            project_id (str): project id associated with the species
            human_orthologs_go_ontology_url: ontology for go sentences for human orthologs, to be used in case of
                information poor genes
            human_orthologs_go_associations_url: association file url for go sentences for human orthologs, to be used
                in case of information poor genes
            expression_cluster_anatomy_prefix (str): prefix for files on WB ftp site related to expression cluster
                anatomy
            expression_cluster_molreg_prefix (str): prefix for files on WB ftp site related to molecule regulation
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
        self.orthology_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + project_id + \
                             '/annotation/' + species + '.' + project_id + '.' + release_version + '.orthologs.txt.gz'
        self.orthology_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                 project_id, "annotation", species + '.' + project_id + '.' +
                                                 release_version + ".orthologs.txt.gz")
        self.orthologs = defaultdict(lambda: defaultdict(list))
        self.sister_sp_fullname = sister_sp_fullname
        self.protein_domain_url = raw_files_source + '/' + release_version + '/species/' + species + '/' + \
                                  project_id + '/annotation/' + species + '.' + project_id + '.' + release_version + \
                                  '.protein_domains.csv.gz'
        self.protein_domain_cache_path = os.path.join(cache_location, "wormbase", release_version, "species", species,
                                                      project_id, "annotation", species + '.' + project_id +
                                                      '.' + release_version + ".protein_domains.csv.gz")
        self.protein_domains = defaultdict(list)
        self.expression_ontology_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                           "anatomy_ontology." + release_version + ".obo")
        self.expression_ontology_url = raw_files_source + '/' + release_version + '/ONTOLOGY/anatomy_ontology.' + \
                                       release_version + '.obo'
        self.expression_associations_cache_path = os.path.join(cache_location, "wormbase", release_version, "ONTOLOGY",
                                                           "anatomy_association." + release_version + ".wb")
        self.expression_associations_url = raw_files_source + '/' + release_version + \
                                           '/ONTOLOGY/anatomy_association.' + release_version + '.wb'
        if expression_cluster_anatomy_prefix:
            self.expression_cluster_anatomy_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                                  release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                                  "/" + expression_cluster_anatomy_prefix + \
                                                  "ECsummary_anatomy." + release_version + ".txt"
            self.expression_cluster_anatomy_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                                      "ExprClusterSummary",
                                                                      expression_cluster_anatomy_prefix +
                                                                      "ECsummary_anatomy." + release_version + ".txt")
            self.expression_cluster_anatomy_data = defaultdict(list)
        else:
            self.expression_cluster_anatomy_data = None
        if expression_cluster_molreg_prefix:
            self.expression_cluster_molreg_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                                 release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                                 "/" + expression_cluster_molreg_prefix + "ECsummary_molReg." + \
                                                 release_version + ".txt"
            self.expression_cluster_molreg_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                                     "ExprClusterSummary",
                                                                     expression_cluster_molreg_prefix +
                                                                     "ECsummary_molReg." + release_version + ".txt")
            self.expression_cluster_molreg_data = defaultdict(list)
        else:
            self.expression_cluster_molreg_data = None

    def load_gene_data_from_file(self) -> None:
        """load gene list from pre-set file location"""
        if not self.gene_data or len(self.gene_data.items()) == 0:
            self.gene_data = {}
            file_path = self._get_cached_file(cache_path=self.gene_data_cache_path, file_source_url=self.gene_data_url)
            with open(file_path) as file:
                for line in file:
                    fields = line.strip().split(',')
                    name = fields[2] if fields[2] != '' else fields[3]
                    self.gene_data["WB:" + fields[1]] = Gene("WB:" + fields[1], name, fields[4] == "Dead", False)

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, config: GenedescConfigParser) -> None:
        if associations_type == DataType.GO or associations_type == DataType.EXPR:
            super().load_associations_from_file(associations_type=associations_type, associations_url=associations_url,
                                                associations_cache_path=associations_cache_path, config=config)
            if associations_type == DataType.EXPR:
                associations = []
                primary_ids = set()
                for subj_associations in self.expression_associations.associations_by_subj.values():
                    for association in subj_associations:
                        if len(association["qualifiers"]) == 0 or "Partial" in association["qualifiers"] or "Certain" \
                                in association["qualifiers"]:
                            association["qualifiers"] = ["Verified"]
                            associations.append(association)
                            primary_ids.add(association["object"]["id"])
                        elif association["object"]["id"] not in primary_ids:
                            associations.append(association)
                self.expression_associations = AssociationSetFactory().create_from_assocs(
                    assocs=associations, ontology=self.expression_ontology)
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
                            gene_id = linearr[2]
                            if linearr[1] == "allele":
                                gene_id = linearr[4].split(",")[0]
                            associations.append({"source_line": line,
                                                 "subject": {
                                                     "id": gene_id,
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
                                                                       terms_blacklist=config.get_module_property(
                                                                           module=Module.DO_EXP_AND_BIO,
                                                                           prop=ConfigModuleProperty.EXCLUDE_TERMS))

    def load_orthology_from_file(self):
        orthology_file = self._get_cached_file(cache_path=self.orthology_cache_path,
                                               file_source_url=self.orthology_url)
        orthologs = defaultdict(list)
        gene_id = ""
        header = True
        for line in open(orthology_file):
            if not line.startswith("#"):
                if line.strip() == "=":
                    header = True
                    self.orthologs["WB:" + gene_id] = orthologs
                    orthologs = defaultdict(list)
                elif header:
                    gene_id = line.strip().split()[0]
                    header = False
                else:
                    ortholog_arr = line.strip().split("\t")
                    orthologs[ortholog_arr[0]].append(ortholog_arr[1:4])

    def get_best_orthologs_for_gene(self, gene_id: str, orth_species_full_name: List[str],
                                    sister_species_data_fetcher: DataManager = None,
                                    ecode_priority_list: List[str] = None):
        """get the list of best orthologs for a gene

        given a gene and a list of possible species from which the best orthologs for the gene can come from, this
        function returns the list of best orthologs for the first species in the list that have orthology data for the
        gene. If a sister species data manager is provided, the list of orthologs for that species is sorted by number
        of annotations

        Args:
            gene_id (str): the gene id
            orth_species_full_name (List[str]): list of species from which to look for the best orthologs
            sister_species_data_fetcher: data fetcher for
            ecode_priority_list:

        Returns (Tuple[List, str]): a tuple containing the list of best orthologs and the name of the organism from
            which they are taken

        """
        best_orthologs = None
        curr_orth_fullname = None
        if len(orth_species_full_name) > 0:
            for curr_orth_fullname in orth_species_full_name:
                if curr_orth_fullname in self.orthologs[gene_id]:
                    orthologs = self.orthologs[gene_id][curr_orth_fullname]
                    # for human orthologs, take only those predicted by more than 1 method
                    if len(orth_species_full_name) == 1 and orth_species_full_name[0] == "Homo sapiens":
                        orthologs = [ortholog for ortholog in orthologs if len(ortholog[2].split(";")) > 1]
                    orthologs_keys = []
                    if len(orthologs) > 0:
                        if len(orthologs) > 1:
                            for ortholog in orthologs:
                                if sister_species_data_fetcher:
                                    orthologs_keys.append([ortholog[0], ortholog[1], len(ortholog[2].split(";")),
                                                           len(sister_species_data_fetcher.get_annotations_for_gene(
                                                               gene_id=ortholog[0], annot_type=DataType.GO,
                                                               priority_list=ecode_priority_list))])
                                else:
                                    orthologs_keys.append([ortholog[0], ortholog[1], len(ortholog[2].split(";"))])
                            if sister_species_data_fetcher:
                                best_orthologs = [sorted(orthologs_keys, key=lambda x: (x[2], x[3]), reverse=True)[0][0:2]]
                            else:
                                best_orthologs = [[orth_key[0], orth_key[1]] for orth_key in
                                                  sorted(orthologs_keys, key=lambda x: x[2], reverse=True) if
                                                  orth_key[2] == max([orth[2] for orth in orthologs_keys])]
                        else:
                            best_orthologs = [[orthologs[0][0], orthologs[0][1]]]
                        break
        return best_orthologs, curr_orth_fullname

    def load_protein_domain_information(self):
        """load protein domain data"""
        protein_domain_file = self._get_cached_file(cache_path=self.protein_domain_cache_path,
                                                    file_source_url=self.protein_domain_url)
        for line in open(protein_domain_file):
            linearr = line.strip().split("\t")
            if len(linearr) > 3 and linearr[3] != "":
                self.protein_domains[linearr[0]] = [domain[0:-1].split(" \"") if len(domain[0:-1].split(" \"")) > 1 else
                                                    [domain, ""] for domain in linearr[3:]]

    def _load_expression_cluster_file(self, file_cache_path, file_url, load_into_data):
        expr_clust_file = self._get_cached_file(cache_path=file_cache_path, file_source_url=file_url)
        header = True
        for line in open(expr_clust_file):
            if not header:
                linearr = line.strip().split("\t")
                load_into_data[linearr[0]] = linearr[1:]
                load_into_data[linearr[0]][2] = self.transform_expression_cluster_terms(load_into_data[
                                                                                            linearr[0]][2].split(","))
                if load_into_data[linearr[0]] and load_into_data[linearr[0]][3]:
                    load_into_data[linearr[0]][3] = [word for study in load_into_data[linearr[0]][3].split(",") for
                                                     word in study.split(" ") if word != "study" and word != "analysis"]
            else:
                header = False

    def load_expression_cluster_data(self):
        """load all expression cluster data"""
        if self.expression_cluster_anatomy_data is not None:
            self._load_expression_cluster_file(self.expression_cluster_anatomy_cache_path,
                                               self.expression_cluster_anatomy_url,
                                               self.expression_cluster_anatomy_data)
        if self.expression_cluster_molreg_data is not None:
            self._load_expression_cluster_file(self.expression_cluster_molreg_cache_path,
                                               self.expression_cluster_molreg_url,
                                               self.expression_cluster_molreg_data)

    def get_expression_cluster_feature(self, gene_id, expression_cluster_type: ExpressionClusterType,
                                       feature: ExpressionClusterFeature):
        """get an expression cluster feature for a specific type

        Args:
            gene_id (str): the id of the gene
            expression_cluster_type (ExpressionClusterType): the type of expression data to read
            feature (ExpressionClusterFeature): the feature to read

        Returns: the requested feature
        """
        target = None
        if expression_cluster_type == ExpressionClusterType.ANATOMY:
            target = self.expression_cluster_anatomy_data
        elif expression_cluster_type == ExpressionClusterType.MOLREG:
            target = self.expression_cluster_molreg_data
        idx = 2
        if feature == ExpressionClusterFeature.TERMS:
            idx = 2
        elif feature == ExpressionClusterFeature.STUDIES:
            idx = 3
        if target and gene_id in target and target[gene_id][idx]:
            return target[gene_id][idx]
        return None

    @staticmethod
    def transform_expression_cluster_terms(terms_list: List[str]):
        inflect_engine = inflect.engine()
        return ["the " + term if inflect_engine.singular_noun(term.split(" ")[-1]) is False else
                term for term in terms_list]

    def load_all_data_from_file(self, config: GenedescConfigParser) -> None:
        """load all data types from pre-set file locations

        Args:
            config (GenedescConfigParser): configuration object where to read properties
        """
        self.load_gene_data_from_file()
        self.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=self.go_ontology_url,
                                     ontology_cache_path=self.go_ontology_cache_path,
                                     config=config)
        self.load_associations_from_file(associations_type=DataType.GO, associations_url=self.go_associations_url,
                                         associations_cache_path=self.go_associations_cache_path, config=config)
        self.load_ontology_from_file(ontology_type=DataType.DO, ontology_url=self.do_ontology_url,
                                     ontology_cache_path=self.do_ontology_cache_path, config=config)
        self.load_associations_from_file(associations_type=DataType.DO, associations_url=self.do_associations_url,
                                         associations_cache_path=self.do_associations_cache_path,
                                         config=config)
        self.load_orthology_from_file()
        self.load_protein_domain_information()
