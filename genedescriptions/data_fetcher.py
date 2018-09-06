import gzip
import urllib.request
import shutil
import os
import re
import inflect

from collections import defaultdict
from typing import List, Iterable, Dict
from ontobio import AssociationSetFactory
from genedescriptions.descriptions_rules import set_all_depths_in_subgraph, Gene, DataType, \
    is_human_ortholog_name_valid, rename_human_ortholog_name
from ontobio.ontol_factory import OntologyFactory
from ontobio.ontol import Ontology
from ontobio.assocmodel import AssociationSet
from ontobio.io.gafparser import GafParser


class DataFetcher(object):
    """retrieve data for gene descriptions from different sources"""

    def __init__(self, go_relations: List[str] = None, do_relations: List[str] = None, use_cache: bool = False):
        """create a new a data fetcher

        Args:
            go_relations (List[str]): list of ontology relations to be used for GO
            do_relations (List[str]): list of ontology relations to be used for DO
            use_cache (bool): whether to use cached files
        """
        self.go_associations = None
        self.go_ontology = None
        self.do_ontology = None
        self.do_associations = None
        self.gene_data = None
        self.expression_ontology = None
        self.expression_associations = None
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

        Args:
            include_dead_genes (bool): whether to include dead genes in the results
            include_pseudo_genes (bool): whether to include pseudo genes in the results
        Returns:
            Gene: data for one gene per each call, including gene_id and gene_name
        """
        if self.gene_data and len(self.gene_data) > 0:
            for gene_id, gene_obj in self.gene_data.items():
                if (include_dead_genes or not gene_obj.dead) and (include_pseudo_genes or not gene_obj.pseudo):
                    yield gene_obj

    @staticmethod
    def remove_blacklisted_annotations(association_set: AssociationSet, ontology: Ontology,
                                       terms_blacklist: List[str] = None) -> AssociationSet:
        """remove annotations linked to blacklisted ontology terms from an association set

        Args:
            association_set (AssociationSet): the original association set
            ontology (Ontology): the ontology linked to the annotations
            terms_blacklist (List[str]): the list of ontology terms related to the annotations to be removed
        Returns:
            AssociationSet: the filtered annotations
        """
        if terms_blacklist:
            associations = []
            for subj_associations in association_set.associations_by_subj.values():
                for association in subj_associations:
                    if association["object"]["id"] not in terms_blacklist:
                        associations.append(association)
            return AssociationSetFactory().create_from_assocs(assocs=associations, ontology=ontology)
        else:
            return association_set

    @staticmethod
    def rename_ontology_terms(ontology: Ontology, terms_replacement_regex: Dict[str, str] = None) -> None:
        """rename ontology terms based on regular expression matching

        Args:
            ontology (Ontology): the ontology containing the terms to be renamed
            terms_replacement_regex (Dict[str, str]): a dictionary containing the regular expression to be applied for
                renaming terms. Each key must be a regular expression to search for terms and the associated value
                another regular expression that defines the final result
        """
        if terms_replacement_regex:
            for regex_to_substitute, regex_target in terms_replacement_regex.items():
                for node in ontology.search(regex_to_substitute, is_regex=True):
                    ontology.node(node)["label"] = re.sub(regex_to_substitute, regex_target,
                                                          ontology.node(node)["label"])

    def set_ontology(self, ontology_type: DataType, ontology: Ontology,
                     terms_replacement_regex: Dict[str, str] = None) -> None:
        """set the go ontology and apply terms renaming

        Args:
            ontology_type (DataType): the type of ontology to set
            ontology (Ontology): an ontology object to set as go ontology
            terms_replacement_regex (Dict[str, str]): a dictionary containing the regular expression to be applied for
                renaming terms. Each key must be a regular expression to search for terms and the associated value
                another regular expression that defines the final result
        """
        new_ontology = None
        if ontology_type == DataType.GO:
            self.go_ontology = ontology.subontology(relations=self.go_relations)
            new_ontology = self.go_ontology
        elif ontology_type == DataType.DO:
            self.do_ontology = ontology.subontology(relations=self.do_relations)
            new_ontology = self.do_ontology
        elif ontology_type == DataType.EXPR:
            self.expression_ontology = ontology.subontology()
            new_ontology = self.expression_ontology
        self.rename_ontology_terms(ontology=new_ontology, terms_replacement_regex=terms_replacement_regex)
        for root_id in new_ontology.get_roots():
            set_all_depths_in_subgraph(ontology=new_ontology, root_id=root_id, relations=None)

    def load_ontology_from_file(self, ontology_type: DataType, ontology_url: str, ontology_cache_path: str,
                                terms_replacement_regex: Dict[str, str] = None) -> None:
        """load go ontology from file

        Args:
            ontology_type (DataType): the type of ontology to set
            ontology_url (str): url to the ontology file
            ontology_cache_path (str): path to cache file for the ontology
            terms_replacement_regex (Dict[str, str])]: a dictionary containing the regular expression to be applied for
                renaming terms. Each key must be a regular expression to search for terms and the associated value
                another regular expression that defines the final result
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
        elif ontology_type == DataType.EXPR:
            self.expression_ontology = OntologyFactory().create(self._get_cached_file(
                file_source_url=ontology_url, cache_path=ontology_cache_path)).subontology()
            new_ontology = self.expression_ontology
        if terms_replacement_regex:
            self.rename_ontology_terms(ontology=new_ontology, terms_replacement_regex=terms_replacement_regex)
        if ontology_type == DataType.EXPR:
            inflect_engine = inflect.engine()
            for term in self.expression_ontology.nodes():
                if "label" in self.expression_ontology.node(term) and \
                        inflect_engine.singular_noun(self.expression_ontology.node(term)["label"].split(" ")[-1]) is \
                        False:
                    self.expression_ontology.node(term)["label"] = "the " + self.expression_ontology.node(term)["label"]
        for root_id in new_ontology.get_roots():
            set_all_depths_in_subgraph(ontology=new_ontology, root_id=root_id, relations=None)

    def set_associations(self, associations_type: DataType, associations: AssociationSet,
                         exclusion_list: List[str] = None) -> None:
        """set the go annotations and remove blacklisted annotations

        Args:
            associations_type (DataType): the type of associations to set
            associations (AssociationSet): an association object to set as go annotations
            exclusion_list (List[str]): the list of ontology terms related to the annotations to be removed
        """
        if associations_type == DataType.GO:
            self.go_associations = self.remove_blacklisted_annotations(association_set=associations,
                                                                       ontology=self.go_ontology,
                                                                       terms_blacklist=exclusion_list)
        elif associations_type == DataType.DO:
            self.do_associations = self.remove_blacklisted_annotations(association_set=associations,
                                                                       ontology=self.do_ontology,
                                                                       terms_blacklist=exclusion_list)
        elif associations_type == DataType.EXPR:
            self.expression_associations = self.remove_blacklisted_annotations(association_set=associations,
                                                                               ontology=self.do_ontology,
                                                                               terms_blacklist=exclusion_list)

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, exclusion_list: List[str]) -> None:
        """load go associations from file

        Args:
            associations_type (DataType): the type of associations to set
            associations_url (str): url to the association file
            associations_cache_path (str): path to cache file for the associations
            exclusion_list (List[str]): the list of ontology terms related to the annotations to be removed
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
        elif associations_type == DataType.EXPR:
            self.expression_associations = AssociationSetFactory().create_from_assocs(assocs=GafParser().parse(
                file=self._get_cached_file(cache_path=associations_cache_path, file_source_url=associations_url),
                skipheader=True), ontology=self.expression_ontology)
            self.expression_associations = self.remove_blacklisted_annotations(
                association_set=self.expression_associations, ontology=self.expression_ontology,
                terms_blacklist=exclusion_list)

    def get_annotations_for_gene(self, gene_id: str, annot_type: DataType = DataType.GO,
                                 include_obsolete: bool = False, include_negative_results: bool = False,
                                 priority_list: Iterable = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "IC", "ISS",
                                                            "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD",
                                                            "RCA", "IEA")) -> List[Dict]:
        """
        retrieve go annotations for a given gene id and a given type. The annotations are unique for each pair
        <gene_id, term_id>. This means that when multiple annotations for the same pair are found in the go data, the
        one with the evidence code with highest priority is returned (see the *priority_list* parameter to set the
        priority according to evidence codes)

        Args:
            gene_id (str): the id of the gene related to the annotations to retrieve, in standard format
            annot_type (DataType): type of annotations to read
            include_obsolete (bool): whether to include obsolete annotations
            include_negative_results (bool): whether to include negative results
            priority_list (List[str]): the priority list for the evidence codes. If multiple annotations with the same
                term are found, only the one with highest priority is returned. The first element in the list has the
                highest priority, whereas the last has the lowest. Only annotations with evidence codes in the priority
                list are returned. All other annotations are ignored
                of annotations for the gene
        Returns:
            List[Dict]: the list of annotations for the given gene
        """
        dataset = None
        ontology = None
        if annot_type == DataType.GO:
            dataset = self.go_associations
            ontology = self.go_ontology
        elif annot_type == DataType.DO:
            dataset = self.do_associations
            ontology = self.do_ontology
        elif annot_type == DataType.EXPR:
            dataset = self.expression_associations
            ontology = self.expression_ontology
        if dataset is not None and ontology is not None:
            priority_map = dict(zip(priority_list, reversed(range(len(list(priority_list))))))
            annotations = [annotation for annotation in dataset.associations(gene_id) if (include_obsolete or
                                                                                          not ontology.is_obsolete(
                                                                                              annotation["object"]["id"]))
                           and (include_negative_results or "NOT" not in annotation["qualifiers"])]
            id_selected_annotation = {}
            for annotation in annotations:
                if annotation["evidence"]["type"] in priority_map.keys():
                    if annotation["object"]["id"] in id_selected_annotation:
                        if priority_map[annotation["evidence"]["type"]] > \
                                priority_map[id_selected_annotation[annotation["object"]["id"]]["evidence"]["type"]]:
                            id_selected_annotation[annotation["object"]["id"]] = annotation
                    else:
                        id_selected_annotation[annotation["object"]["id"]] = annotation
            return [annotation for annotation in id_selected_annotation.values()]
        else:
            return []

    def set_gene_data(self, gene_data: List[Gene]):
        for gene in gene_data:
            self.gene_data[gene.id] = gene

    def load_gene_data_from_file(self):
        pass

    @staticmethod
    def get_human_gene_props(use_ensembl_id: bool = True):
        """ retrieve data for human genes, including Ensembl ID, symbol, name, and family symbol and name

        Args:
            use_ensembl_id (bool): use ensembl id as key instead of hgnc id
        Returns:
            Dict[str, List[str]]: a dictionary of all human genes properties, indexed by Ensembl ID

        """
        human_genes_props = defaultdict(list)
        human_content_w_ensmbl = urllib.request.urlopen("https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col="
                                                        "gd_pub_ensembl_id&status=Approved&status=Entry+Withdrawn&statu"
                                                        "s_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgn"
                                                        "c_dbtag=on&submit=submit")
        human_content_w_fam_sym = urllib.request.urlopen(
            "https://www.genenames.org/cgi-bin/genefamilies/download-all/tsv")

        header = True
        for line in human_content_w_ensmbl:
            if not header:
                linearr = line.decode("utf-8").split("\t")
                linearr[-1] = linearr[-1].strip()
                if linearr[1] != "":
                    human_genes_props[linearr[0][5:]] = [linearr[1]]
            else:
                header = False
        header = True
        for line in human_content_w_fam_sym:
            if not header:
                linearr = line.decode("utf-8").split("\t")
                linearr[-1] = linearr[-1].strip()
                if is_human_ortholog_name_valid(linearr[2]):
                    human_genes_props[linearr[0]].extend([linearr[1], rename_human_ortholog_name(linearr[2]),
                                                          linearr[9], linearr[10]])
                else:
                    del human_genes_props[linearr[0]]
            else:
                header = False
        if use_ensembl_id:
            return {v[0]: v[1:] for k, v in human_genes_props.items()}
        else:
            return {k: v[1:] for k, v in human_genes_props.items()}

    @staticmethod
    def get_ensembl_hgnc_ids_map():
        human_genes_props = {}
        human_content_w_ensmbl = urllib.request.urlopen("https://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col="
                                                        "gd_pub_ensembl_id&status=Approved&status=Entry+Withdrawn&statu"
                                                        "s_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgn"
                                                        "c_dbtag=on&submit=submit")
        header = True
        for line in human_content_w_ensmbl:
            if not header:
                linearr = line.decode("utf-8").split("\t")
                linearr[-1] = linearr[-1].strip()
                if linearr[1] != "":
                    human_genes_props[linearr[1]] = linearr[0]
            else:
                header = False
        return human_genes_props


class WBDataFetcher(DataFetcher):
    """data fetcher for WormBase raw files for a single species"""

    def __init__(self, raw_files_source: str, cache_location: str, release_version: str, species: str, project_id: str,
                 go_relations: List[str] = None, do_relations: List[str] = None, use_cache: bool = False,
                 sister_sp_fullname: str = "", human_orthologs_go_ontology_url: str = None,
                 human_orthologs_go_associations_url: str = None):
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
        self.expression_enriched_extra_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                             release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                             "/ceECsummary_anatomy." + release_version + ".txt"
        self.expression_enriched_extra_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                                 "ExprClusterSummary", "ceECsummary_anatomy." +
                                                                 release_version + ".txt")
        self.expression_enriched_extra_data = defaultdict(list)
        self.expression_enriched_bma_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                            release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                            "/bmaECsummary_anatomy." + release_version + ".txt"
        self.expression_enriched_bma_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                               "ExprClusterSummary", "bmaECsummary_anatomy." +
                                                               release_version + ".txt")
        self.expression_enriched_bma_data = defaultdict(list)
        self.expression_affected_bma_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                            release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                            "/bmaECsummary_molReg." + release_version + ".txt"
        self.expression_affected_bma_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                               "ExprClusterSummary", "bmaECsummary_molReg." +
                                                               release_version + ".txt")
        self.expression_affected_bma_data = defaultdict(list)
        self.expression_enriched_ppa_url = "ftp://caltech.wormbase.org/pub/wormbase/ExprClusterSummary/" + \
                                           release_version[0:-1] + str(int(release_version[-1]) + 1) + \
                                           "/ppaECsummary_anatomy." + release_version + ".txt"
        self.expression_enriched_ppa_cache_path = os.path.join(cache_location, "wormbase", release_version,
                                                               "ExprClusterSummary", "ppaECsummary_anatomy." +
                                                               release_version + ".txt")
        self.expression_enriched_ppa_data = defaultdict(list)

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
                                    associations_cache_path: str, exclusion_list: List[str]) -> None:
        if associations_type == DataType.GO or associations_type == DataType.EXPR:
            super().load_associations_from_file(associations_type=associations_type, associations_url=associations_url,
                                                associations_cache_path=associations_cache_path,
                                                exclusion_list=exclusion_list)
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
                                                                       terms_blacklist=exclusion_list)

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
                                    sister_species_data_fetcher: DataFetcher = None,
                                    ecode_priority_list: List[str] = None):
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

    def get_all_orthologs_for_gene(self, gene_id: str, organism: str):
        return self.orthologs[gene_id][organism]

    def load_protein_domain_information(self):
        protein_domain_file = self._get_cached_file(cache_path=self.protein_domain_cache_path,
                                                    file_source_url=self.protein_domain_url)
        for line in open(protein_domain_file):
            linearr = line.strip().split("\t")
            if len(linearr) > 3 and linearr[3] != "":
                self.protein_domains[linearr[0]] = [domain[0:-1].split(" \"") if len(domain[0:-1].split(" \"")) > 1 else
                                                    [domain, ""] for domain in linearr[3:]]

    def load_expression_enriched_extra_info(self):
        expr_enriched_extra_file = self._get_cached_file(cache_path=self.expression_enriched_extra_cache_path,
                                                         file_source_url=self.expression_enriched_extra_url)
        header = True
        for line in open(expr_enriched_extra_file):
            if not header:
                linearr = line.strip().split("\t")
                self.expression_enriched_extra_data[linearr[0]] = [word for study in linearr[4].split(",") for word in
                                                                   study.split(" ") if word != "study" and word !=
                                                                   "analysis"]
            else:
                header = False

    def load_bma_expression_data(self):
        expr_enriched_extra_file = self._get_cached_file(cache_path=self.expression_enriched_bma_cache_path,
                                                         file_source_url=self.expression_enriched_bma_url)
        header = True
        for line in open(expr_enriched_extra_file):
            if not header:
                linearr = line.strip().split("\t")
                self.expression_enriched_bma_data[linearr[0]] = linearr[1:]
                self.expression_enriched_bma_data[linearr[0]][2] = \
                    self.transform_expression_cluster_terms(self.expression_enriched_bma_data[linearr[0]][2].split(","))
                if self.expression_enriched_bma_data[linearr[0]] and self.expression_enriched_bma_data[linearr[0]][3]:
                    self.expression_enriched_bma_data[linearr[0]][3] = \
                        [word for study in self.expression_enriched_bma_data[linearr[0]][3].split(",") for word in
                         study.split(" ") if word != "study" and word != "analysis"]
            else:
                header = False
        expr_enriched_extra_file = self._get_cached_file(cache_path=self.expression_affected_bma_cache_path,
                                                         file_source_url=self.expression_affected_bma_url)
        header = True
        for line in open(expr_enriched_extra_file):
            if not header:
                linearr = line.strip().split("\t")
                self.expression_affected_bma_data[linearr[0]] = linearr[1:]
                self.expression_affected_bma_data[linearr[0]][2] = self.transform_expression_cluster_terms(
                    self.expression_affected_bma_data[linearr[0]][2].split(","))
                if self.expression_affected_bma_data[linearr[0]] and self.expression_affected_bma_data[linearr[0]][3]:
                    self.expression_affected_bma_data[linearr[0]][3] = \
                        [word for study in self.expression_affected_bma_data[linearr[0]][3].split(",") for word in
                         study.split(" ") if word != "study"]
            else:
                header = False

    def load_ppa_expression_data(self):
        expr_enriched_extra_file = self._get_cached_file(cache_path=self.expression_enriched_ppa_cache_path,
                                                         file_source_url=self.expression_enriched_ppa_url)
        header = True
        for line in open(expr_enriched_extra_file):
            if not header:
                linearr = line.strip().split("\t")
                self.expression_enriched_ppa_data[linearr[0]] = linearr[1:]
                self.expression_enriched_ppa_data[linearr[0]][2] = self.transform_expression_cluster_terms(
                    self.expression_enriched_ppa_data[linearr[0]][2].split(","))
                if self.expression_enriched_ppa_data[linearr[0]][3]:
                    self.expression_enriched_ppa_data[linearr[0]][3] = \
                        [word for study in self.expression_enriched_ppa_data[linearr[0]][3].split(",") for word in
                         study.split(" ") if word != "study" and word != "analysis"]
            else:
                header = False

    @staticmethod
    def transform_expression_cluster_terms(terms_list: List[str]):
        inflect_engine = inflect.engine()
        return ["the " + term if inflect_engine.singular_noun(term.split(" ")[-1]) is False else
                term for term in terms_list]

    def load_all_data_from_file(self, go_terms_replacement_regex: Dict[str, str] = None,
                                go_terms_exclusion_list: List[str] = None,
                                do_terms_replacement_regex: Dict[str, str] = None,
                                do_terms_exclusion_list: List[str] = None) -> None:
        """load all data types from pre-set file locations

        Args:
            go_terms_replacement_regex (Dict[str, str]): a dictionary containing the regular expression to be applied
                for renaming terms. Each key must be a regular expression to search for terms and the associated value
                another regular expression that defines the final result
            go_terms_exclusion_list (List[str]): the list of GO ontology terms related to the GO annotations to be
                removed
            do_terms_replacement_regex (Dict[str, str]): a dictionary containing the regular expression to be applied
                for renaming DO terms. Each key must be a regular expression to search for terms and the associated
                value another regular expression that defines the final result
            do_terms_exclusion_list (List[str]): the list of DO ontology terms related to the DO annotations to be
                removed
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
        self.load_orthology_from_file()
        self.load_protein_domain_information()

