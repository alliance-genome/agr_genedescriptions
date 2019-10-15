import logging
import os
import re

import inflect

from collections import defaultdict
from typing import List
from ontobio import AssociationSetFactory
from ontobio.io.assocparser import AssocParserConfig
from ontobio.io.gafparser import GafParser
from genedescriptions.commons import DataType, Gene, Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import ExpressionClusterFeature, DataManager, ExpressionClusterType


logger = logging.getLogger("WB Data Manager")


class WBDataManager(DataManager):
    """data fetcher for WormBase raw files for a single species"""

    def __init__(self, config: GenedescConfigParser, species: str, go_relations: List[str] = None,
                 do_relations: List[str] = None, use_cache: bool = False):
        """create a new data fetcher for WormBase. Files will be downloaded from WB ftp site. For convenience, file
        locations are automatically generated and stored in class variables ending in _url for remote filed and
        _cache_path for caching

        Args:
            species (str): WormBase species to fetch
        """
        self.config = config
        raw_files_source = config.get_wb_raw_file_sources()
        cache_location = config.get_cache_dir()
        release_version = config.get_wb_release()
        organisms_info = config.get_wb_organisms_info()
        project_id = organisms_info[species]["project_id"]
        self.sister_sp_fullname = ""
        if "main_sister_species" in organisms_info[species] and "full_name" in \
                organisms_info[organisms_info[species]["main_sister_species"]]:
            self.sister_sp_fullname = organisms_info[organisms_info[species]["main_sister_species"]]["full_name"]
        self.orth_fullnames = ""
        if "ortholog" in organisms_info[species] and all(["full_name" in organisms_info[ortholog_sp] for ortholog_sp in
                                                          organisms_info[species]["ortholog"]]):
            self.orth_fullnames = [organisms_info[ortholog_sp]["full_name"] for ortholog_sp in
                                   organisms_info[species]["ortholog"]]
        expression_cluster_anatomy_prefix = organisms_info[species]["ec_anatomy_prefix"] if \
            "ec_anatomy_prefix" in organisms_info[species] else None
        expression_cluster_molreg_prefix = organisms_info[species]["ec_molreg_prefix"] if \
            "ec_molreg_prefix" in organisms_info[species] else None
        expression_cluster_genereg_prefix = organisms_info[species]["ec_genereg_prefix"] if \
            "ec_genereg_prefix" in organisms_info[species] else None
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
        self.expression_cluster_anatomy_url = self._get_expression_cluster_url(
            prefix=expression_cluster_anatomy_prefix, ec_type="anatomy", release_version=release_version)
        self.expression_cluster_anatomy_cache_path = self._get_expression_cluster_cache_path(
            prefix=expression_cluster_anatomy_prefix, ec_type="anatomy", release_version=release_version,
            cache_location=cache_location)
        self.expression_cluster_anatomy_data = defaultdict(list) if self.expression_cluster_anatomy_url else None
        self.expression_cluster_molreg_url = self._get_expression_cluster_url(
            prefix=expression_cluster_molreg_prefix, ec_type="molReg", release_version=release_version)
        self.expression_cluster_molreg_cache_path = self._get_expression_cluster_cache_path(
            prefix=expression_cluster_molreg_prefix, ec_type="molReg", release_version=release_version,
            cache_location=cache_location)
        self.expression_cluster_molreg_data = defaultdict(list) if self.expression_cluster_molreg_url else None
        self.expression_cluster_genereg_url = self._get_expression_cluster_url(
            prefix=expression_cluster_genereg_prefix, ec_type="geneReg", release_version=release_version)
        self.expression_cluster_genereg_cache_path = self._get_expression_cluster_cache_path(
            prefix=expression_cluster_genereg_prefix, ec_type="geneReg", release_version=release_version,
            cache_location=cache_location)
        self.expression_cluster_genereg_data = defaultdict(list) if self.expression_cluster_genereg_url else None

    @staticmethod
    def _get_expression_cluster_url(prefix, ec_type, release_version):
        if prefix:
            return "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + release_version[0:2] + \
                   str(int(release_version[2:]) + 1) + "/" + prefix + "ECsummary_" + ec_type + "." + release_version + \
                   ".txt"
        else:
            return None

    @staticmethod
    def _get_expression_cluster_cache_path(prefix, ec_type, release_version, cache_location):
        if prefix:
            return os.path.join(cache_location, "wormbase", release_version, "ExprClusterSummary", prefix +
                                "ECsummary_" + ec_type + "." + release_version + ".txt")
        else:
            return None

    def load_gene_data_from_file(self) -> None:
        """load gene list from pre-set file location"""
        logger.info("Loading genes data from file")
        if not self.gene_data or len(self.gene_data.items()) == 0:
            self.gene_data = {}
            file_path = self._get_cached_file(cache_path=self.gene_data_cache_path, file_source_url=self.gene_data_url)
            with open(file_path) as file:
                for line in file:
                    fields = line.strip().split(',')
                    if fields[1].startswith("WBGene"):
                        name = fields[2] if fields[2] != '' else fields[3]
                        self.gene_data["WB:" + fields[1]] = Gene("WB:" + fields[1], name, fields[4] == "Dead", False)

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, config: GenedescConfigParser,
                                    association_additional_url: str = None,
                                    association_additional_cache_path: str = None) -> None:
        logger.info("Loading associations from file")
        assoc_config = AssocParserConfig(remove_double_prefixes=True, paint=True)
        if associations_type == DataType.GO:
            super().load_associations_from_file(associations_type=associations_type, associations_url=associations_url,
                                                associations_cache_path=associations_cache_path, config=config)
        elif associations_type == DataType.EXPR:
            associations = []
            file_path = self._get_cached_file(cache_path=associations_cache_path, file_source_url=associations_url)
            for line in open(file_path):
                if not line.strip().startswith("!"):
                    linearr = line.strip().split("\t")
                    if self.expression_ontology.node(linearr[4]):
                        gene_id = linearr[0] + ":" + linearr[1]
                        qualifiers = linearr[3].split("|")
                        if len(qualifiers) == 0 or "Partial" in qualifiers or "Certain" in qualifiers:
                            qualifiers = ["Verified"]
                        associations.append(DataManager.create_annotation_record(
                            line, gene_id, linearr[2], linearr[11], linearr[12], linearr[4], qualifiers, linearr[8],
                            linearr[6], linearr[5].split("|"), linearr[14], linearr[13]))
            self.expression_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                                      ontology=self.expression_ontology)
            self.expression_associations = self.remove_blacklisted_annotations(
                association_set=self.expression_associations, ontology=self.expression_ontology,
                terms_blacklist=config.get_module_property(module=Module.EXPRESSION,
                                                           prop=ConfigModuleProperty.EXCLUDE_TERMS))
        elif associations_type == DataType.DO:
            self.do_associations = AssociationSetFactory().create_from_assocs(
                assocs=GafParser(config=assoc_config).parse(file=self._get_cached_file(
                    cache_path=associations_cache_path, file_source_url=associations_url), skipheader=True),
                ontology=self.do_ontology)
            if association_additional_cache_path and association_additional_url:
                associations = []
                for subj_associations in self.do_associations.associations_by_subj.values():
                    for association in subj_associations:
                        if association["evidence"]["type"] == "IEA":
                            associations.append(association)
                file_path = self._get_cached_file(cache_path=association_additional_cache_path,
                                                  file_source_url=association_additional_url)
                header = True
                for line in open(file_path):
                    if not line.strip().startswith("!"):
                        if not header:
                            linearr = line.strip().split("\t")
                            if self.do_ontology.node(linearr[10]) and "IEA" not in linearr[16] and \
                                    linearr[1] in ["gene", "allele"]:
                                gene_ids = ["WB:" + linearr[3]]
                                if linearr[1] == "allele":
                                    gene_ids = linearr[4].split(",")
                                if linearr[8] == "is_marker_for":
                                    linearr[16] = "BMK"
                                elif linearr[8] == "is_implicated_in" or linearr[8] == "is_model_of":
                                    linearr[16] = "IMP"
                                for gene_id in gene_ids:
                                    associations.append(DataManager.create_annotation_record(
                                        line, gene_id, linearr[3], linearr[1], linearr[0], linearr[10],
                                        linearr[9].split("|"), "D", linearr[16], linearr[18].split("|"),
                                        linearr[20], linearr[19]))
                        else:
                            header = False
                self.do_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                                  ontology=self.do_ontology)
            self.do_associations = self.remove_blacklisted_annotations(
                association_set=self.do_associations, ontology=self.do_ontology,
                terms_blacklist=config.get_module_property(module=Module.DO_EXPERIMENTAL,
                                                           prop=ConfigModuleProperty.EXCLUDE_TERMS))

    def load_orthology_from_file(self):
        logger.info("Loading orthology from file")
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
                    if not ortholog_arr[1].startswith("PRJEB28388"):
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
        logger.info("Getting list of best orthologs for gene")
        best_orthologs = None
        curr_orth_fullname = None
        if len(orth_species_full_name) > 0:
            for curr_orth_fullname in orth_species_full_name:
                if curr_orth_fullname in self.orthologs[gene_id]:
                    orthologs = self.orthologs[gene_id][curr_orth_fullname]
                    # for human orthologs, take only those predicted by more than 1 method - removed
                    # if len(orth_species_full_name) == 1 and orth_species_full_name[0] == "Homo sapiens":
                    #     orthologs = [ortholog for ortholog in orthologs if len(ortholog[2].split(";")) > 1]
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
                                best_orthologs = [sorted(orthologs_keys, key=lambda x: (x[2], x[3]),
                                                         reverse=True)[0][0:2]]
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
        logger.info("Loading protein domain information from file")
        protein_domain_file = self._get_cached_file(cache_path=self.protein_domain_cache_path,
                                                    file_source_url=self.protein_domain_url)
        for line in open(protein_domain_file):
            linearr = line.strip().split("\t")
            if len(linearr) > 3 and linearr[3] != "":
                self.protein_domains[linearr[0]] = [domain[0:-1].split(" \"") if len(domain[0:-1].split(" \"")) > 1 else
                                                    [domain, ""] for domain in linearr[3:]]

    @staticmethod
    def get_replaced_terms_arr(terms, terms_replacement_regex):
        new_terms = terms
        for regex_to_substitute, regex_target in terms_replacement_regex.items():
            tmp_terms = []
            for term in new_terms:
                tmp_terms.append(re.sub(regex_to_substitute, regex_target, term))
            new_terms = tmp_terms
        return new_terms

    def _load_expression_cluster_file(self, file_cache_path, file_url, load_into_data,
                                      add_to_expression_ontology_annotations: bool = False):
        expr_clust_file = self._get_cached_file(cache_path=file_cache_path, file_source_url=file_url)
        header = True
        associations = []
        terms_ids_map = {}
        if add_to_expression_ontology_annotations:
            associations = [association for subj_associations in
                            self.expression_associations.associations_by_subj.values() for association in
                            subj_associations]
        terms_replacement_regex = self.config.get_module_property(module=Module.EXPRESSION,
                                                                  prop=ConfigModuleProperty.RENAME_TERMS)
        for line in open(expr_clust_file):
            if not header:
                linearr = line.strip().split("\t")
                load_into_data[linearr[0]] = linearr[1:]
                load_into_data[linearr[0]][2] = WBDataManager.get_replaced_terms_arr(
                    load_into_data[linearr[0]][2].split(","), terms_replacement_regex)
                if load_into_data[linearr[0]] and load_into_data[linearr[0]][3]:
                    load_into_data[linearr[0]][3] = [word.replace(" study", "").replace(" analysis", "") for word in
                                                     load_into_data[linearr[0]][3].split(",")]
                if add_to_expression_ontology_annotations:
                    for term in load_into_data[linearr[0]][2]:
                        if term not in terms_ids_map:
                            term_ids = self.expression_ontology.resolve_names([term])
                            if term_ids:
                                terms_ids_map[term] = term_ids[0]
                            else:
                                terms_ids_map[term] = None
                        if term in terms_ids_map and terms_ids_map[term]:
                            associations.append(DataManager.create_annotation_record(
                                line, "WB:" + linearr[0], "", "gene", "", terms_ids_map[term], ["Enriched"], "A", "IDA",
                                "", "", ""))
            else:
                header = False
        if add_to_expression_ontology_annotations:
            self.set_associations(DataType.EXPR, associations=AssociationSetFactory().create_from_assocs(
                assocs=associations, ontology=self.expression_ontology), config=self.config)

    def load_expression_cluster_data(self):
        """load all expression cluster data"""
        logger.info("Loading expression cluster data from file")
        if self.expression_cluster_anatomy_data is not None:
            self._load_expression_cluster_file(self.expression_cluster_anatomy_cache_path,
                                               self.expression_cluster_anatomy_url,
                                               self.expression_cluster_anatomy_data,
                                               add_to_expression_ontology_annotations=True)
        if self.expression_cluster_molreg_data is not None:
            self._load_expression_cluster_file(self.expression_cluster_molreg_cache_path,
                                               self.expression_cluster_molreg_url,
                                               self.expression_cluster_molreg_data)
        if self.expression_cluster_genereg_data is not None:
            self._load_expression_cluster_file(self.expression_cluster_genereg_cache_path,
                                               self.expression_cluster_genereg_url,
                                               self.expression_cluster_genereg_data)

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
        elif expression_cluster_type == ExpressionClusterType.GENEREG:
            target = self.expression_cluster_genereg_data
        idx = 2
        if feature == ExpressionClusterFeature.TERMS:
            idx = 2
        elif feature == ExpressionClusterFeature.STUDIES:
            idx = 3
        if target and gene_id in target and target[gene_id][idx]:
            return target[gene_id][idx]
        return None

    def load_all_data_from_file(self) -> None:
        """load all data types from pre-set file locations"""
        self.load_gene_data_from_file()
        self.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=self.go_ontology_url,
                                     ontology_cache_path=self.go_ontology_cache_path,
                                     config=self.config)
        self.load_associations_from_file(associations_type=DataType.GO, associations_url=self.go_associations_url,
                                         associations_cache_path=self.go_associations_cache_path, config=self.config)
        self.load_ontology_from_file(ontology_type=DataType.DO, ontology_url=self.do_ontology_url,
                                     ontology_cache_path=self.do_ontology_cache_path, config=self.config)
        self.load_associations_from_file(associations_type=DataType.DO, associations_url=self.do_associations_url,
                                         associations_cache_path=self.do_associations_cache_path,
                                         association_additional_cache_path=self.do_associations_new_cache_path,
                                         association_additional_url=self.do_associations_new_url, config=self.config)
        self.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url=self.expression_ontology_url,
                                     ontology_cache_path=self.expression_ontology_cache_path, config=self.config)
        self.load_associations_from_file(associations_type=DataType.EXPR,
                                         associations_url=self.expression_associations_url,
                                         associations_cache_path=self.expression_associations_cache_path,
                                         config=self.config)
        self.load_orthology_from_file()
        self.load_expression_cluster_data()
        self.load_protein_domain_information()
