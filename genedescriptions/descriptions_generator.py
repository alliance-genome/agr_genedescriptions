from collections import ChainMap
from typing import Set

import inflect
import re

from genedescriptions.commons import Sentence, Module, DataType, TrimmingResult
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from genedescriptions.ontology_tools import *
from genedescriptions.sentence_generation_functions import _get_single_sentence, compose_sentence
from genedescriptions.trimming import CONF_TO_TRIMMING_CLASS_DEFAULT
from genedescriptions.trimming_custom import CONF_TO_TRIMMING_CLASS_CUSTOM

logger = logging.getLogger(__name__)


CONF_TO_TRIMMING_CLASS = ChainMap(CONF_TO_TRIMMING_CLASS_DEFAULT, CONF_TO_TRIMMING_CLASS_CUSTOM)


class ModuleSentences(object):
    def __init__(self, sentences):
        self.sentences = sentences

    def get_description(self):
        return " and ".join([sentence.text for sentence in self.sentences])

    def get_ids(self, experimental_only: bool = False):
        return [term_id for sentence in self.sentences for term_id in sentence.terms_ids if not experimental_only or
                sentence.evidence_group.startswith("EXPERIMENTAL")]

    def contains_sentences(self):
        return len(self.sentences) > 0


class SentenceMerger(object):
    def __init__(self):
        self.postfix_list = []
        self.terms_ids = set()
        self.term_postfix_dict = {}
        self.evidence_groups = []
        self.term_evgroup_dict = {}
        self.additional_prefix = ""
        self.aspect = ""
        self.qualifier = ""
        self.ancestors_covering_multiple_terms = set()
        self.any_trimmed = False


class OntologySentenceGenerator(object):
    """generates sentences based on description rules"""

    def __init__(self, gene_id: str, module: Module, data_manager: DataManager, config: GenedescConfigParser,
                 limit_to_group: str = None, humans: bool = False):
        """initialize sentence generator object

        Args:
            config (GenedescConfigParser): an optional config object from which to read the options
            limit_to_group (str): limit the evidence codes to the specified group
        """
        annot_type = None
        if module == Module.DO_ORTHOLOGY or module == Module.DO_EXPERIMENTAL or module == module.DO_BIOMARKER:
            self.ontology = data_manager.do_ontology
            annot_type = DataType.DO
        elif module == Module.GO:
            self.ontology = data_manager.go_ontology
            annot_type = DataType.GO
        elif module == Module.EXPRESSION:
            self.ontology = data_manager.expression_ontology
            annot_type = DataType.EXPR
        self.evidence_groups_priority_list = config.get_evidence_groups_priority_list(module=module)
        self.prepostfix_sentences_map = config.get_prepostfix_sentence_map(module=module, humans=humans)
        self.terms_groups = defaultdict(lambda: defaultdict(set))
        ev_codes_groups_maps = config.get_evidence_codes_groups_map(module=module)
        annotations = data_manager.get_annotations_for_gene(gene_id=gene_id, annot_type=annot_type,
                                                            priority_list=config.get_annotations_priority(
                                                                module=module))
        self.annotations = annotations
        self.module = module
        self.data_manager = data_manager
        self.annot_type = annot_type
        evidence_codes_groups_map = {evcode: group for evcode, group in ev_codes_groups_maps.items() if
                                     limit_to_group is None or limit_to_group in ev_codes_groups_maps[evcode]}
        prepostfix_special_cases_sent_map = config.get_prepostfix_sentence_map(module=module, special_cases_only=True,
                                                                               humans=humans)
        self.cat_several_words = config.get_module_property(module=self.module,
                                                            prop=ConfigModuleProperty.CUTOFF_SEVERAL_CATEGORY_WORD)
        if not self.cat_several_words:
            self.cat_several_words = {'F': 'functions', 'P': 'processes', 'C': 'components', 'D': 'diseases',
                                      'A': 'tissues'}
        self.del_overlap = config.get_module_property(module=self.module, prop=ConfigModuleProperty.REMOVE_OVERLAP)
        self.remove_parents = config.get_module_property(module=self.module,
                                                         prop=ConfigModuleProperty.DEL_PARENTS_IF_CHILD)
        self.remove_child_terms = config.get_module_property(module=self.module,
                                                             prop=ConfigModuleProperty.DEL_CHILDREN_IF_PARENT)
        self.max_terms = config.get_module_property(module=self.module,
                                                    prop=ConfigModuleProperty.MAX_NUM_TERMS_IN_SENTENCE)
        self.exclude_terms = config.get_module_property(module=self.module, prop=ConfigModuleProperty.EXCLUDE_TERMS)
        if not self.exclude_terms:
            self.exclude_terms = []
        self.cutoff_final_word = config.get_module_property(module=self.module,
                                                            prop=ConfigModuleProperty.CUTOFF_SEVERAL_WORD)
        self.rename_cell = config.get_module_property(module=self.module, prop=ConfigModuleProperty.RENAME_CELL)
        self.terms_already_covered = set()
        self.dist_root = config.get_module_property(module=self.module,
                                                    prop=ConfigModuleProperty.DISTANCE_FROM_ROOT)
        if not self.dist_root:
            self.dist_root = {'F': 1, 'P': 1, 'C': 2, 'D': 3, 'A': 3}
        self.add_mul_common_anc = config.get_module_property(
            module=self.module, prop=ConfigModuleProperty.ADD_MULTIPLE_TO_COMMON_ANCEST)
        self.trimming_algorithm = config.get_module_property(module=self.module,
                                                             prop=ConfigModuleProperty.TRIMMING_ALGORITHM)
        self.slim_set = self.data_manager.get_slim(module=self.module)
        self.slim_bonus_perc = config.get_module_property(module=self.module, prop=ConfigModuleProperty.SLIM_BONUS_PERC)
        self.add_mul_common_anc = config.get_module_property(module=self.module,
                                                             prop=ConfigModuleProperty.ADD_MULTIPLE_TO_COMMON_ANCEST)
        self.nodeids_blacklist = config.get_module_property(module=self.module, prop=ConfigModuleProperty.EXCLUDE_TERMS)
        self.config = config
        if len(annotations) > 0:
            for annotation in annotations:
                if annotation["evidence"]["type"] in evidence_codes_groups_map:
                    aspect = annotation["aspect"]
                    ev_group = evidence_codes_groups_map[annotation["evidence"]["type"]]
                    qualifier = "_".join(sorted(annotation["qualifiers"])) if "qualifiers" in annotation else ""
                    if prepostfix_special_cases_sent_map and (aspect, ev_group, qualifier) in \
                            prepostfix_special_cases_sent_map:
                        for special_case in prepostfix_special_cases_sent_map[(aspect, ev_group, qualifier)]:
                            if re.match(re.escape(special_case[1]), self.ontology.label(annotation["object"]["id"],
                                                                                        id_if_null=True)):
                                ev_group = evidence_codes_groups_map[annotation["evidence"]["type"]] + \
                                           str(special_case[0])
                                if ev_group not in self.evidence_groups_priority_list:
                                    self.evidence_groups_priority_list.insert(self.evidence_groups_priority_list.index(
                                        evidence_codes_groups_map[annotation["evidence"]["type"]]) + 1, ev_group)
                                break
                    self.terms_groups[(aspect, qualifier)][ev_group].add(annotation["object"]["id"])

    def get_module_sentences(self, aspect: str, qualifier: str = '',
                             keep_only_best_group: bool = False, merge_groups_with_same_prefix: bool = False,
                             high_priority_term_ids: List[str] = None):
        """generate description for a specific combination of aspect and qualifier

        Args:
            aspect (str): a data type aspect
            qualifier (str): qualifier
            keep_only_best_group (bool): whether to get only the evidence group with highest priority and discard
                the other evidence groups
            merge_groups_with_same_prefix (bool): whether to merge the phrases for evidence groups with the same prefix
            high_priority_term_ids (List[str]): list of ids for terms that must always appear in the sentence with
                higher priority than the other terms. Trimming is not applied to these terms
        Returns:
            ModuleSentences: the module sentences
        """
        sentences = []
        evidence_group_priority = {eg: p for p, eg in enumerate(self.evidence_groups_priority_list)}
        for terms, evidence_group, priority in sorted([(t, eg, evidence_group_priority[eg]) for eg, t in
                                                       self.terms_groups[(aspect, qualifier)].items()],
                                                      key=lambda x: x[2]):
            trimming_result = self.reduce_num_terms(terms=terms, min_distance_from_root=self.dist_root[aspect],
                                                    high_priority_term_ids=high_priority_term_ids)
            if (aspect, evidence_group, qualifier) in self.prepostfix_sentences_map and len(terms) > 0:
                sentences.append(
                    _get_single_sentence(
                        node_ids=trimming_result.final_terms, ontology=self.ontology, aspect=aspect,
                        evidence_group=evidence_group, qualifier=qualifier,
                        prepostfix_sentences_map=self.prepostfix_sentences_map,
                        terms_merged=False, trimmed=trimming_result.trimming_applied,
                        add_others=trimming_result.partial_coverage,
                        truncate_others_generic_word=self.cutoff_final_word,
                        truncate_others_aspect_words=self.cat_several_words,
                        ancestors_with_multiple_children=trimming_result.multicovering_nodes,
                        rename_cell=self.rename_cell, config=self.config,
                        put_anatomy_male_at_end=True if aspect == 'A' else False))
                if keep_only_best_group:
                    return ModuleSentences(sentences)
        if merge_groups_with_same_prefix:
            sentences = self.merge_sentences_with_same_prefix(
                sentences=sentences, remove_parent_terms=self.remove_parents, rename_cell=self.rename_cell,
                high_priority_term_ids=high_priority_term_ids, put_anatomy_male_at_end=True if aspect == 'A' else False)
        return ModuleSentences(sentences)

    def reduce_num_terms(self, terms: Set[str], min_distance_from_root: int = 0,
                         high_priority_term_ids: List[str] = None) -> TrimmingResult:
        """
        Reduce the initial set of terms by resolving parent child relationships, deleting overlap with previous
        sections, and trimming, if needed, the final set of terms

        Args:
            terms (Set[str]): initial set of terms to reduce
            min_distance_from_root (int): minimum distance from root to use for trimming
            high_priority_term_ids (List[str]): list of terms with high priority to be treated separately

        Returns:
            TrimmingResult: the reduced set of terms with additional information on the nature of the terms
        """
        trimming_result = TrimmingResult()
        if self.del_overlap:
            terms -= self.terms_already_covered
        if self.exclude_terms:
            terms -= set(self.exclude_terms)
        if self.remove_parents:
            terms = OntologySentenceGenerator.remove_parents_if_child_present(
                terms, self.ontology, self.terms_already_covered, high_priority_term_ids)
        if 0 < self.max_terms < len(terms):
            trimming_result = self.trim_with_priority(terms, min_distance_from_root, high_priority_term_ids)
        else:
            trimming_result.final_terms = terms
            self.terms_already_covered.update(terms)
        if self.remove_child_terms:
            trimming_result.final_terms = self.remove_children_if_parents_present(
                terms=trimming_result.final_terms, ontology=self.ontology,
                terms_already_covered=self.terms_already_covered, high_priority_terms=high_priority_term_ids,
                ancestors_covering_multiple_children=trimming_result.multicovering_nodes if
                self.add_mul_common_anc else None)
        return trimming_result

    def trim_with_priority(self, terms: Set[str], min_distance_from_root: int = 0,
                           high_priority_terms: List[str] = None) -> TrimmingResult:
        """
        Trim terms based on common ancestors and keep terms with high priority in a separate set, to be trimmed
        separately

        Args:
            terms (Set[str]): the set of terms to trim
            min_distance_from_root (int): minimum distance from root for selecting common ancestors
            high_priority_terms (List[str]): list of term IDs with higher priority

        Returns:
            TrimmingResult: the trimmed terms, with their additional information
        """
        comb_trim_res = TrimmingResult()
        trim_res_hp = TrimmingResult()
        trim_res_lp = TrimmingResult()
        add_multicover_nodes = set()
        trim_res_hp.final_terms = [t for t in terms if t in high_priority_terms] if high_priority_terms else []
        if len(trim_res_hp.final_terms) > self.max_terms:
            trim_res_hp.final_terms = self.remove_children_if_parents_present(
                terms=trim_res_hp.final_terms, ontology=self.ontology, terms_already_covered=self.terms_already_covered,
                ancestors_covering_multiple_children=add_multicover_nodes if self.add_mul_common_anc else None)
        if len(trim_res_hp.final_terms) > self.max_terms:
            logger.debug("Reached maximum number of terms. Applying trimming to high priority terms")
            trim_res_hp = self.trim_terms(trim_res_hp.final_terms, min_distance_from_root)
        trim_res_lp.final_terms = [t for t in terms if t not in high_priority_terms] if high_priority_terms else terms
        trimming_threshold = self.max_terms - len(trim_res_hp.final_terms)
        if 0 < trimming_threshold < len(trim_res_lp.final_terms):
            trim_res_lp = self.trim_terms(trim_res_lp.final_terms, min_distance_from_root)
        comb_trim_res.partial_coverage = trim_res_hp.partial_coverage or trim_res_lp.partial_coverage or \
                                         trimming_threshold <= 0 < len(trim_res_lp.final_terms)
        comb_trim_res.final_terms = trim_res_hp.final_terms
        # remove exact overlap
        comb_trim_res.final_terms.extend(list(set(trim_res_lp.final_terms) - set(trim_res_hp.final_terms)))
        # cutoff terms - if number of terms with high priority is higher than max_num_terms
        comb_trim_res.final_terms = comb_trim_res.final_terms[0:self.max_terms]
        if self.add_mul_common_anc:
            comb_trim_res.multicovering_nodes = trim_res_hp.multicovering_nodes | trim_res_lp.multicovering_nodes | \
                                                add_multicover_nodes
        comb_trim_res.trimming_applied = trim_res_hp.trimming_applied or trim_res_lp.trimming_applied
        comb_trim_res.covered_nodes = trim_res_hp.covered_nodes | trim_res_lp.covered_nodes
        return comb_trim_res

    def trim_terms(self, terms: List[str], min_dist_from_root: int = 0) -> TrimmingResult:
        return CONF_TO_TRIMMING_CLASS[self.trimming_algorithm](
            ontology=self.ontology, min_distance_from_root=min_dist_from_root, nodeids_blacklist=self.nodeids_blacklist,
            slim_terms_ic_bonus_perc=self.slim_bonus_perc, slim_set=self.slim_set).process(node_ids=list(terms),
                                                                                           max_num_nodes=self.max_terms)

    @staticmethod
    def remove_children_if_parents_present(terms, ontology, terms_already_covered: Set[str] = None,
                                           high_priority_terms: List[str] = None,
                                           ancestors_covering_multiple_children: Set[str] = None):
        terms_nochildren = []
        for term in terms:
            if len(set(ontology.ancestors(term)).intersection(set(terms))) == 0 or (high_priority_terms and
                                                                                    term in high_priority_terms):
                terms_nochildren.append(term)
            elif ancestors_covering_multiple_children is not None:
                ancestors_covering_multiple_children.update({ontology.label(term_id, id_if_null=True) for term_id in
                                                             set(ontology.ancestors(term)).intersection(set(terms))})
        if len(terms_nochildren) < len(terms):
            if terms_already_covered is not None:
                terms_already_covered.update(set(terms) - set(terms_nochildren))
            logger.debug("Removed " + str(len(terms) - len(terms_nochildren)) + " children from terms")
            return terms_nochildren
        else:
            return terms

    @staticmethod
    def remove_parents_if_child_present(terms, ontology, terms_already_covered: Set[str] = None,
                                        high_priority_terms: List[str] = None):
        terms_no_ancestors = list(set(terms) - set([ancestor for node_id in terms for ancestor in
                                                    ontology.ancestors(node_id) if not high_priority_terms or
                                                    ancestor not in high_priority_terms]))
        if len(terms) > len(terms_no_ancestors):
            if terms_already_covered is not None:
                terms_already_covered.update(set(terms) - set(terms_no_ancestors))
            logger.debug("Removed " + str(len(terms) - len(terms_no_ancestors)) + " parents from terms")
            return terms_no_ancestors
        else:
            return terms

    def merge_sentences_with_same_prefix(self, sentences: List[Sentence], remove_parent_terms: bool = True,
                                         rename_cell: bool = False, high_priority_term_ids: List[str] = None,
                                         put_anatomy_male_at_end: bool = False):
        """merge sentences with the same prefix

        Args:
            sentences (List[Sentence]): a list of sentences
            remove_parent_terms (bool): whether to remove parent terms if present in the merged set of terms
            rename_cell (bool): whether to rename the term 'cell'
            high_priority_term_ids (List[str]): list of ids for terms that must always appear in the sentence with
                higher priority than the other terms. Trimming is not applied to these terms
        Returns:
            List[Sentence]: the list of merged sentences, sorted by (merged) evidence group priority
        """
        merged_sentences = defaultdict(SentenceMerger)
        for sentence in sentences:
            prefix = self.prepostfix_sentences_map[(sentence.aspect, sentence.evidence_group, sentence.qualifier)][0]
            merged_sentences[prefix].postfix_list.append(self.prepostfix_sentences_map[(sentence.aspect,
                                                                                        sentence.evidence_group,
                                                                                        sentence.qualifier)][1])
            merged_sentences[prefix].aspect = sentence.aspect
            merged_sentences[prefix].qualifier = sentence.qualifier
            merged_sentences[prefix].terms_ids.update(sentence.terms_ids)
            for term in sentence.terms_ids:
                merged_sentences[prefix].term_postfix_dict[term] = self.prepostfix_sentences_map[
                    (sentence.aspect, sentence.evidence_group, sentence.qualifier)][1]
            merged_sentences[prefix].evidence_groups.append(sentence.evidence_group)
            for term in sentence.terms_ids:
                merged_sentences[prefix].term_evgroup_dict[term] = sentence.evidence_group
            if sentence.additional_prefix:
                merged_sentences[prefix].additional_prefix = sentence.additional_prefix
            merged_sentences[prefix].ancestors_covering_multiple_terms.update(
                sentence.ancestors_covering_multiple_terms)
            merged_sentences[prefix].any_trimmed = merged_sentences[prefix].any_trimmed or sentence.trimmed
        if remove_parent_terms:
            for prefix, sent_merger in merged_sentences.items():
                terms_no_ancestors = sent_merger.terms_ids - set([ancestor for node_id in sent_merger.terms_ids for
                                                                  ancestor in self.ontology.ancestors(node_id) if not
                                                                  high_priority_term_ids or ancestor not in
                                                                  high_priority_term_ids])
                if len(sent_merger.terms_ids) > len(terms_no_ancestors):
                    logger.debug("Removed " + str(len(sent_merger.terms_ids) - len(terms_no_ancestors)) +
                                 " parents from terms while merging sentences with same prefix")
                    sent_merger.terms_ids = terms_no_ancestors
        return [Sentence(prefix=prefix, terms_ids=list(sent_merger.terms_ids),
                         postfix=OntologySentenceGenerator.merge_postfix_phrases(sent_merger.postfix_list),
                         text=compose_sentence(prefix=prefix,
                                               term_names=[self.ontology.label(node, id_if_null=True) for node in
                                                           sent_merger.terms_ids],
                                               postfix=OntologySentenceGenerator.merge_postfix_phrases(
                                                   sent_merger.postfix_list),
                                               additional_prefix=sent_merger.additional_prefix,
                                               ancestors_with_multiple_children=sent_merger.ancestors_covering_multiple_terms,
                                               rename_cell=rename_cell, config=self.config,
                                               put_anatomy_male_at_end=put_anatomy_male_at_end),
                         aspect=sent_merger.aspect, evidence_group=", ".join(sent_merger.evidence_groups),
                         terms_merged=True, trimmed=sent_merger.any_trimmed,
                         additional_prefix=sent_merger.additional_prefix, qualifier=sent_merger.qualifier,
                         ancestors_covering_multiple_terms=sent_merger.ancestors_covering_multiple_terms)
                for prefix, sent_merger in merged_sentences.items() if len(sent_merger.terms_ids) > 0]

    @staticmethod
    def merge_postfix_phrases(postfix_phrases: List[str]) -> str:
        """merge postfix phrases and remove possible redundant text at the beginning at at the end of the phrases

        Args:
            postfix_phrases (List[str]): the phrases to merge
        Returns:
            str: the merged postfix phrase
        """
        postfix_phrases = [postfix for postfix in postfix_phrases if postfix]
        if postfix_phrases and len(postfix_phrases) > 0:
            if len(postfix_phrases) > 1:
                inf_engine = inflect.engine()
                shortest_phrase = sorted(zip(postfix_phrases, [len(phrase) for phrase in postfix_phrases]),
                                         key=lambda x: x[1])[0][0]
                first_part = ""
                for idx, letter in enumerate(shortest_phrase):
                    if all(map(lambda x: x[idx] == shortest_phrase[idx], postfix_phrases)):
                        first_part += letter
                    else:
                        break
                last_part = ""
                for idx, letter in zip(range(len(shortest_phrase)), reversed([l for l in shortest_phrase])):
                    if all(map(lambda x: x[len(x) - idx - 1] == shortest_phrase[len(shortest_phrase) - idx - 1],
                               postfix_phrases)):
                        last_part = letter + last_part
                    else:
                        break
                new_phrases = [phrase.replace(first_part, "").replace(last_part, "") for phrase in postfix_phrases]
                if len(last_part.strip().split(" ")) == 1:
                    last_part = inf_engine.plural(last_part)
                if len(new_phrases) > 2:
                    return first_part + ", ".join(new_phrases[0:-1]) + ", and " + new_phrases[-1] + last_part
                elif len(new_phrases) > 1:
                    return first_part + " and ".join(new_phrases) + last_part
                else:
                    return first_part + new_phrases[0] + last_part
            else:
                return postfix_phrases[0]
        else:
            return ""


