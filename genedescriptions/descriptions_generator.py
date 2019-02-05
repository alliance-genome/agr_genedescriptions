import inflect
import re

from genedescriptions.commons import Sentence, Module, DataType
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager
from genedescriptions.ontology_tools import *
from genedescriptions.sentence_generation_functions import _get_single_sentence, compose_sentence


logger = logging.getLogger(__name__)


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

    def get_module_sentences(self,  config: GenedescConfigParser, aspect: str, qualifier: str = '',
                             keep_only_best_group: bool = False, merge_groups_with_same_prefix: bool = False,
                             high_priority_term_ids: List[str] = None):
        """generate description for a specific combination of aspect and qualifier

        Args:
            config (GenedescConfigParser): a configuration object from which to read properties
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
        cat_several_words = config.get_module_property(module=self.module,
                                                       prop=ConfigModuleProperty.CUTOFF_SEVERAL_CATEGORY_WORD)
        del_overlap = config.get_module_property(module=self.module, prop=ConfigModuleProperty.REMOVE_OVERLAP)
        remove_parents = config.get_module_property(module=self.module,
                                                    prop=ConfigModuleProperty.DEL_PARENTS_IF_CHILD)
        remove_child_terms = config.get_module_property(module=self.module,
                                                        prop=ConfigModuleProperty.DEL_CHILDREN_IF_PARENT)
        max_terms = config.get_module_property(module=self.module,
                                               prop=ConfigModuleProperty.MAX_NUM_TERMS_IN_SENTENCE)
        exclude_terms = config.get_module_property(module=self.module, prop=ConfigModuleProperty.EXCLUDE_TERMS)
        cutoff_final_word = config.get_module_property(module=self.module,
                                                       prop=ConfigModuleProperty.CUTOFF_SEVERAL_WORD)
        rename_cell = config.get_module_property(module=self.module, prop=ConfigModuleProperty.RENAME_CELL)
        if not cat_several_words:
            cat_several_words = {'F': 'functions', 'P': 'processes', 'C': 'components', 'D': 'diseases', 'A': 'tissues'}
        sentences = []
        terms_already_covered = set()
        evidence_group_priority = {eg: p for p, eg in enumerate(self.evidence_groups_priority_list)}
        for terms, evidence_group, priority in sorted([(t, eg, evidence_group_priority[eg]) for eg, t in
                                                       self.terms_groups[(aspect, qualifier)].items()],
                                                      key=lambda x: x[2]):
            ancestors_covering_multiple_children = set()
            if del_overlap:
                terms -= terms_already_covered
            if exclude_terms:
                terms -= set(exclude_terms)
            if remove_parents:
                terms_no_ancestors = terms - set([ancestor for node_id in terms for ancestor in
                                                  self.ontology.ancestors(node_id) if not high_priority_term_ids or
                                                  ancestor not in high_priority_term_ids])
                if len(terms) > len(terms_no_ancestors):
                    logger.debug("Removed " + str(len(terms) - len(terms_no_ancestors)) + " parents from terms")
                    terms = terms_no_ancestors
            trimmed = False
            add_others = False
            if 0 < max_terms < len(terms):
                trimmed = True
                terms, add_others, ancestors_covering_multiple_children = self.get_trimmed_terms_by_common_ancestor(
                    terms, terms_already_covered, aspect, config, high_priority_term_ids)
            else:
                terms_already_covered.update(terms)
            if remove_child_terms:
                terms = [term for term in terms if
                         len(set(self.ontology.ancestors(term)).intersection(set(terms))) == 0 or (
                                 high_priority_term_ids and term in high_priority_term_ids)]
            if (aspect, evidence_group, qualifier) in self.prepostfix_sentences_map and len(terms) > 0:
                sentences.append(
                    _get_single_sentence(
                        node_ids=terms, ontology=self.ontology, aspect=aspect, evidence_group=evidence_group,
                        qualifier=qualifier, prepostfix_sentences_map=self.prepostfix_sentences_map,
                        terms_merged=False, trimmed=trimmed, add_others=add_others,
                        truncate_others_generic_word=cutoff_final_word,
                        truncate_others_aspect_words=cat_several_words,
                        ancestors_with_multiple_children=ancestors_covering_multiple_children, rename_cell=rename_cell))
                if keep_only_best_group:
                    return ModuleSentences(sentences)
        if merge_groups_with_same_prefix:
            sentences = self.merge_sentences_with_same_prefix(
                sentences=sentences, remove_parent_terms=remove_parents, rename_cell=rename_cell,
                high_priority_term_ids=high_priority_term_ids)
        return ModuleSentences(sentences)

    def get_trimmed_terms_by_common_ancestor(self, terms: Set[str], terms_already_covered, aspect: str,
                                             config: GenedescConfigParser, high_priority_terms: List[str] = None):
        dist_root = config.get_module_property(module=self.module, prop=ConfigModuleProperty.DISTANCE_FROM_ROOT)
        add_mul_common_anc = config.get_module_property(module=self.module,
                                                        prop=ConfigModuleProperty.ADD_MULTIPLE_TO_COMMON_ANCEST)
        max_terms = config.get_module_property(module=self.module,
                                               prop=ConfigModuleProperty.MAX_NUM_TERMS_IN_SENTENCE)
        trimming_algorithm = config.get_module_property(module=self.module,
                                                        prop=ConfigModuleProperty.TRIMMING_ALGORITHM)
        slim_set = self.data_manager.get_slim(module=self.module)
        slim_bonus_perc = config.get_module_property(module=self.module, prop=ConfigModuleProperty.SLIM_BONUS_PERC)
        add_others = False
        ancestors_covering_multiple_children = set()
        if not dist_root:
            dist_root = {'F': 1, 'P': 1, 'C': 2, 'D': 3, 'A': 3}
        terms_high_priority = [term for term in terms if high_priority_terms and term in high_priority_terms]
        if terms_high_priority is None:
            terms_high_priority = []
        terms_already_covered.update(terms_high_priority)
        terms_low_priority = [term for term in terms if not high_priority_terms or term not in
                              high_priority_terms]
        trimming_threshold = max_terms - len(terms_high_priority)
        if 0 < trimming_threshold < len(terms_low_priority):
            merged_terms_coverset = None
            if trimming_algorithm == "naive":
                add_others, merged_terms_coverset = get_best_nodes_naive(
                    node_ids=list(terms_low_priority), ontology=self.ontology, min_distance_from_root=dist_root[aspect])
            elif trimming_algorithm == "ic":
                add_others, merged_terms_coverset = get_best_nodes_ic(
                    node_ids=list(terms_low_priority), ontology=self.ontology, max_number_of_terms=trimming_threshold,
                    min_distance_from_root=dist_root[aspect], slim_terms_ic_bonus_perc=slim_bonus_perc,
                    slim_set=slim_set)
            elif trimming_algorithm == "naive2":
                add_others, merged_terms_coverset = get_best_nodes_lca(
                    node_ids=list(terms_low_priority), ontology=self.ontology, min_distance_from_root=dist_root[aspect])
            if add_mul_common_anc:
                ancestors_covering_multiple_children = {self.ontology.label(term_id, id_if_null=True) for
                                                        term_id, covered_nodes in merged_terms_coverset if
                                                        term_id not in terms_low_priority}
            terms_low_priority = [term_id for term_id, covered_nodes in merged_terms_coverset]
            terms_already_covered.update([e for term_id, covered_nodes in merged_terms_coverset for e in covered_nodes])
        terms = terms_high_priority
        if len(terms) > max_terms:
            # remove children if parent is present in terms for key diseases when they are too many
            terms = [term for term in terms if len(set(self.ontology.ancestors(term)).intersection(set(terms))) == 0]
        if len(terms) > max_terms:
            add_others = True
        terms_low_priority = [term for term in terms_low_priority if term not in terms_high_priority]
        terms.extend(terms_low_priority)
        # cutoff terms - if number of terms with high priority is higher than max_num_terms
        terms = terms[0:max_terms]
        return terms, add_others, ancestors_covering_multiple_children

    def merge_sentences_with_same_prefix(self, sentences: List[Sentence], remove_parent_terms: bool = True,
                                         rename_cell: bool = False, high_priority_term_ids: List[str] = None):
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
                                               rename_cell=rename_cell),
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


