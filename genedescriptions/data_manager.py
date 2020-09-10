import gzip
import logging
import urllib.request
import shutil
import os
import re
import inflect

from enum import Enum
from collections import defaultdict
from typing import List, Iterable, Dict
from ontobio import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory
from ontobio.ontol import Ontology
from ontobio.assocmodel import AssociationSet
from genedescriptions.commons import Gene, DataType, Module, get_module_from_data_type
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.ontology_tools import set_all_depths, set_ic_annot_freq, set_ic_ontology_struct


class ExpressionClusterType(Enum):
    ANATOMY = 1
    MOLREG = 2
    GENEREG = 3


class ExpressionClusterFeature(Enum):
    TERMS = 2
    STUDIES = 3


logger = logging.getLogger(__name__)


class DataManager(object):
    """retrieve data for gene descriptions from different sources"""

    def __init__(self, go_relations: List[str] = None, do_relations: List[str] = None, expr_relations: List[str] = None,
                 use_cache: bool = False):
        """create a new a data fetcher

        Args:
            go_relations (List[str]): list of ontology relations to be used for GO
            do_relations (List[str]): list of ontology relations to be used for DO
            expr_relations (List[str]): list of ontology relations to be used for EXPRESSION
            use_cache (bool): whether to use cached files
        """
        self.go_associations = None
        self.go_ontology = None
        self.do_ontology = None
        self.do_associations = None
        self.gene_data = {}
        self.expression_ontology = None
        self.expression_associations = None
        self.go_relations = go_relations
        self.do_relations = do_relations
        self.expr_relations = expr_relations
        self.go_slim = set()
        self.do_slim = set()
        self.exp_slim = set()
        self.use_cache = use_cache

    def get_ontology(self, data_type: DataType):
        if data_type == DataType.DO:
            return self.do_ontology
        elif data_type == DataType.GO:
            return self.go_ontology
        elif data_type == DataType.EXPR:
            return self.expression_ontology
        return None

    def get_associations(self, data_type: DataType):
        if data_type == DataType.DO:
            return self.do_associations
        elif data_type == DataType.GO:
            return self.go_associations
        elif data_type == DataType.EXPR:
            return self.expression_associations
        return None

    def get_relations(self, data_type: DataType):
        if data_type == DataType.DO:
            return self.do_relations
        elif data_type == DataType.GO:
            return self.go_relations
        elif data_type == DataType.EXPR:
            return self.expr_relations
        return None

    @staticmethod
    def get_slim_cache_path(ontology_cache_path, data_type: DataType):
        slim_name = "slim.obo"
        if data_type == DataType.GO:
            slim_name = "go_slim.obo"
        elif data_type == DataType.DO:
            slim_name = "do_slim.obo"
        elif data_type == DataType.EXPR:
            slim_name = "expr_slim.obo"
        return os.path.join(os.path.dirname(os.path.normpath(ontology_cache_path)), slim_name)

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
        logger.info("Removing blacklisted terms and annotations")
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
        logger.info("Renaming ontology terms")
        if terms_replacement_regex:
            for regex_to_substitute, regex_target in terms_replacement_regex.items():
                for node in ontology.search(regex_to_substitute, is_regex=True):
                    ontology.node(node)["label"] = re.sub(regex_to_substitute, regex_target,
                                                          ontology.node(node)["label"])

    def set_ontology(self, ontology_type: DataType, ontology: Ontology, config: GenedescConfigParser,
                     slim_cache_path: str = None) -> None:
        """set the go ontology and apply terms renaming

        Args:
            ontology_type (DataType): the type of ontology to set
            ontology (Ontology): an ontology object to set as go ontology
            config (GenedescConfigParser): configuration object where to read properties
            slim_cache_path (str): path to slim file to use
        """
        if ontology_type == DataType.GO:
            logger.info("Setting GO ontology")
            if self.go_relations:
                self.go_ontology = ontology.subontology(relations=self.go_relations)
            else:
                self.go_ontology = ontology
        elif ontology_type == DataType.DO:
            logger.info("Setting DO ontology")
            if self.do_relations:
                self.do_ontology = ontology.subontology(relations=self.do_relations)
            else:
                self.do_ontology = ontology
        elif ontology_type == DataType.EXPR:
            logger.info("Setting Expression ontology")
            if self.expr_relations:
                self.expression_ontology = ontology.subontology(relations=self.expr_relations)
            else:
                self.expression_ontology = ontology
        module = get_module_from_data_type(ontology_type)
        ontology = self.get_ontology(data_type=ontology_type)
        terms_replacement_regex = config.get_module_property(module=module, prop=ConfigModuleProperty.RENAME_TERMS)
        if terms_replacement_regex:
            self.rename_ontology_terms(ontology=ontology, terms_replacement_regex=terms_replacement_regex)
        set_all_depths(ontology=ontology, relations=self.get_relations(ontology_type))
        if config.get_module_property(module=module,
                                      prop=ConfigModuleProperty.TRIMMING_ALGORITHM) == "ic":
            set_ic_ontology_struct(ontology=ontology, relations=self.get_relations(ontology_type))
        if slim_cache_path:
            slim_url = config.get_module_property(module=module, prop=ConfigModuleProperty.SLIM_URL)
            self.load_slim(module=module, slim_url=slim_url, slim_cache_path=slim_cache_path)

    @staticmethod
    def add_article_to_nodes(ontology):
        inflect_engine = inflect.engine()
        for term in ontology.nodes():
            if "label" in ontology.node(term) and \
                    inflect_engine.singular_noun(ontology.node(term)["label"].split(" ")[-1]) is False:
                ontology.node(term)["label"] = "the " + ontology.node(term)["label"]

    def load_ontology_from_file(self, ontology_type: DataType, ontology_url: str, ontology_cache_path: str,
                                config: GenedescConfigParser) -> None:
        """load go ontology from file

        Args:
            ontology_type (DataType): the type of ontology to set
            ontology_url (str): url to the ontology file
            ontology_cache_path (str): path to cache file for the ontology
            config (GenedescConfigParser): configuration object where to read properties
        """
        new_ontology = OntologyFactory().create(self._get_cached_file(file_source_url=ontology_url,
                                                                      cache_path=ontology_cache_path))
        self.set_ontology(ontology_type=ontology_type, ontology=new_ontology, config=config,
                          slim_cache_path=self.get_slim_cache_path(ontology_cache_path, ontology_type))

    def load_slim(self, module: Module, slim_url: str, slim_cache_path: str):
        if slim_url and slim_cache_path:
            relations = None
            if module == Module.GO:
                relations = self.go_relations
            elif module == Module.DO_EXPERIMENTAL:
                relations = self.do_relations
            elif module == Module.EXPRESSION:
                relations = None
            slim_onto = OntologyFactory().create(self._get_cached_file(file_source_url=slim_url, cache_path=slim_cache_path)
                                                 ).subontology(relations=relations)
            slim_set = set([node for node in slim_onto.nodes() if "type" in slim_onto.node(node) and
                            slim_onto.node(node)["type"] == "CLASS"])
            if module == Module.GO:
                logger.info("Setting GO Slim")
                self.go_slim = slim_set
            elif module == Module.DO_EXPERIMENTAL:
                logger.info("Setting DO Slim")
                self.do_slim = slim_set
            elif module == Module.EXPRESSION:
                logger.info("Setting Expression Slim")
                self.exp_slim = slim_set

    def get_slim(self, module: Module):
        if module == Module.GO:
            return self.go_slim
        elif module == Module.DO_EXPERIMENTAL:
            return self.do_slim
        elif module == Module.EXPRESSION:
            return self.exp_slim

    @staticmethod
    def remap_associations(associations: AssociationSet, ontology: Ontology, associations_map: Dict[str, str]):
        if not associations_map:
            return associations
        new_associations = []
        for subj_associations in associations.associations_by_subj.values():
            for association in subj_associations:
                if association["object"]["id"] in associations_map:
                    association["object"]["id"] = associations_map[association["object"]["id"]]
                new_associations.append(association)
        return AssociationSetFactory().create_from_assocs(assocs=new_associations, ontology=ontology)

    def set_associations(self, associations_type: DataType, associations: AssociationSet, config: GenedescConfigParser):
        """set the go annotations and remove blacklisted annotations

        Args:
            associations_type (DataType): the type of associations to set
            associations (AssociationSet): an association object to set as go annotations
            config (GenedescConfigParser): configuration object where to read properties
        """
        assocs = self.remap_associations(associations=associations, ontology=self.get_ontology(associations_type),
                                         associations_map=config.get_module_property(
                                             module=get_module_from_data_type(associations_type),
                                             prop=ConfigModuleProperty.REMAP_TERMS))
        assocs = self.remove_blacklisted_annotations(
            association_set=assocs, ontology=self.get_ontology(associations_type),
            terms_blacklist=config.get_module_property(module=get_module_from_data_type(associations_type),
                                                       prop=ConfigModuleProperty.EXCLUDE_TERMS))

        if associations_type == DataType.GO:
            logger.info("Setting GO associations")
            self.go_associations = assocs
        elif associations_type == DataType.DO:
            logger.info("Setting DO associations")
            self.do_associations = assocs
        elif associations_type == DataType.EXPR:
            logger.info("Setting Expression associations")
            self.expression_associations = assocs
        if config.get_module_property(module=get_module_from_data_type(associations_type),
                                      prop=ConfigModuleProperty.TRIMMING_ALGORITHM) == "icGO":
            set_ic_annot_freq(self.get_ontology(associations_type), self.get_associations(associations_type))

    def load_associations_from_file(self, associations_type: DataType, associations_url: str,
                                    associations_cache_path: str, config: GenedescConfigParser) -> None:
        """load go associations from file

        Args:
            associations_type (DataType): the type of associations to set
            associations_url (str): url to the association file
            associations_cache_path (str): path to cache file for the associations
            config (GenedescConfigParser): configuration object where to read properties
        """
        assocs = AssociationSetFactory().create_from_file(file=self._get_cached_file(
            cache_path=associations_cache_path, file_source_url=associations_url),
            ontology=self.get_ontology(associations_type), skim=False)
        self.set_associations(associations_type=associations_type, associations=assocs, config=config)

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
            annotations = [annotation for annotation in dataset.associations(gene_id) if
                           ontology.has_node(annotation["object"]["id"]) and (
                               include_obsolete or ("deprecated" not in ontology.node(annotation["object"]["id"]) or
                                                    not ontology.node(annotation["object"]["id"])["deprecated"]))
                           and (include_negative_results or ("NOT" not in annotation["qualifiers"] and
                                                             not annotation["negated"]))
                           and ontology.label(annotation["object"]["id"])]
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
    def get_human_gene_props():
        """ retrieve data for human genes, including HGNC ID, symbol, and name

        Returns:
            Dict[str, List[str]]: a dictionary of all human genes properties, indexed by HGNC ID

        """
        human_genes_props = defaultdict(list)
        human_content_w_ensmbl = urllib.request.urlopen("https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_pub_ensembl_id&col=gd_app_sym&col=gd_app_name&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit")

        header = True
        for line in human_content_w_ensmbl:
            if not header:
                linearr = line.decode("utf-8").split("\t")
                linearr[-1] = linearr[-1].strip()
                if linearr[1] != "":
                    human_genes_props[linearr[0]] = [linearr[2], linearr[3]]
            else:
                header = False
        return human_genes_props

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

    @staticmethod
    def create_annotation_record(source_line, gene_id, gene_symbol, gene_type, taxon_id, object_id, qualifiers, aspect,
                                 ecode, references, prvdr, date):
        return {"source_line": source_line,
                "subject": {
                    "id": gene_id,
                    "label": gene_symbol,
                    "type": gene_type,
                    "fullname": "",
                    "synonyms": [],
                    "taxon": {"id": taxon_id}
                },
                "object": {
                    "id": object_id,
                    "taxon": ""
                },
                "qualifiers": qualifiers,
                "aspect": aspect,
                "relation": {"id": None},
                "negated": False,
                "evidence": {
                    "type": ecode,
                    "has_supporting_reference": references,
                    "with_support_from": [],
                    "provided_by": prvdr,
                    "date": date
                }}
