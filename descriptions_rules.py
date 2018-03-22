import copy
import logging
import inflect
import re
from collections import namedtuple
from descriptions_writer import SingleDescStats
from ontology_tools import *

GOSentence = namedtuple('GOSentence', ['prefix', 'terms', 'term_ids_dict', 'postfix', 'text', 'go_aspect',
                                       'evidence_group', 'terms_merged', 'additional_prefix'])


class GOSentenceMerger(object):
    def __init__(self):
        self.postfix_list = []
        self.terms = set()
        self.terms_ids_dict = {}
        self.term_postfix_dict = {}
        self.evidence_groups = []
        self.term_evgroup_dict = {}
        self.additional_prefix = ""


class GOSentencesCollection(object):
    """a group of GO sentences indexed by aspect"""
    def __init__(self, evidence_groups_list, go_prepostfix_sentences_map, remove_parent_terms: bool, go_ontology,
                 go_terms_replacement_dict: Dict[str, str] = None):
        self.evidence_groups_list = evidence_groups_list
        self.go_prepostfix_sentences_map = go_prepostfix_sentences_map
        self.sentences_map = {}
        self.remove_parent_terms = remove_parent_terms
        self.go_ontology = go_ontology
        self.go_terms_replacement_dict = go_terms_replacement_dict

    def set_sentence(self, sentence: GOSentence) -> None:
        """add a sentence to the collection

        :param sentence: the sentence to add
        :type sentence: Sentence
        """
        if sentence is not None:
            self.sentences_map[(sentence.go_aspect, sentence.evidence_group)] = sentence

    def get_sentences(self, go_aspect: str, keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False,
                      desc_stats: SingleDescStats = None) -> List[GOSentence]:
        """get all sentences containing the specified aspect

        :param go_aspect: a GO aspect
        :type go_aspect: str
        :param keep_only_best_group: whether to get only the evidence group with highest priority and discard
            the other evidence groups
        :type keep_only_best_group: bool
        :param merge_groups_with_same_prefix: whether to merge the phrases for evidence groups with the same prefix
        :type merge_groups_with_same_prefix: bool
        :param desc_stats: an object containing the description statistics where to save the total number of annotations
            for the gene
        :type desc_stats: SingleDescStats
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[GOSentence]
        """
        sentences = []
        merged_sentences = defaultdict(GOSentenceMerger)
        if desc_stats:
            desc_stats.num_terms_trim_nogroup_priority_nomerge[go_aspect] = 0
            desc_stats.terms_trim_nogroup_priority_nomerge[go_aspect] = []
            for ga, eg in self.sentences_map.keys():
                if ga == go_aspect:
                    desc_stats.num_terms_trim_nogroup_priority_nomerge[go_aspect] += \
                        len(self.sentences_map[(go_aspect, eg)].terms)
                    desc_stats.terms_trim_nogroup_priority_nomerge[go_aspect].extend(
                        self.sentences_map[(go_aspect, eg)].terms)
        terms_in_previous_ev_groups = set()
        for eg in self.evidence_groups_list:
            if (go_aspect, eg) in self.sentences_map:
                if merge_groups_with_same_prefix:
                    prefix = self.go_prepostfix_sentences_map[(go_aspect, eg)][0]
                    merged_sentences[prefix].postfix_list.append(self.go_prepostfix_sentences_map[(go_aspect, eg)][1])
                    merged_sentences[prefix].terms.update([term for term in self.sentences_map[(go_aspect, eg)].terms
                                                          if term not in terms_in_previous_ev_groups])
                    merged_sentences[prefix].terms_ids_dict.update(self.sentences_map[(go_aspect, eg)].term_ids_dict)
                    for term in self.sentences_map[(go_aspect, eg)].terms:
                        merged_sentences[prefix].term_postfix_dict[term] = self.go_prepostfix_sentences_map[
                            (go_aspect, eg)][1]
                    merged_sentences[prefix].evidence_groups.append(eg)
                    for term in self.sentences_map[(go_aspect, eg)].terms:
                        merged_sentences[prefix].term_evgroup_dict[term] = eg
                    terms_in_previous_ev_groups.update(merged_sentences[prefix].terms)
                    if self.sentences_map[(go_aspect, eg)].additional_prefix:
                        merged_sentences[prefix].additional_prefix = \
                            self.sentences_map[(go_aspect, eg)].additional_prefix
                else:
                    sentence_copy = copy.copy(self.sentences_map[(go_aspect, eg)])
                    sentence_copy.terms = [term for term in sentence_copy.terms if term not in
                                           terms_in_previous_ev_groups]
                    sentences.append(sentence_copy)
                    terms_in_previous_ev_groups.update(self.sentences_map[(go_aspect, eg)].terms)
                if keep_only_best_group:
                    break
        if merge_groups_with_same_prefix:
            for prefix, sent_merger in merged_sentences.items():
                # rem parents
                if self.remove_parent_terms:
                    term_ids_no_parents = get_term_ids_without_parents_from_terms_names(
                        go_terms_names=sent_merger.terms, term_ids_dict=sent_merger.terms_ids_dict,
                        ontology=self.go_ontology)
                    if len(sent_merger.terms) > len(term_ids_no_parents):
                        logging.debug("Removed " + str(len(sent_merger.terms) - len(term_ids_no_parents)) +
                                      " parents from terms while merging sentences with same prefix")
                        term_names_no_parents = set([self.go_ontology.query_term(term_id).name for term_id in
                                                     term_ids_no_parents])
                        sent_merger.terms_ids_dict = {key: value for key, value in sent_merger.terms_ids_dict.items()
                                                      if key in set(term_names_no_parents)}
                        sent_merger.terms = list(term_names_no_parents)
                        if self.go_terms_replacement_dict:
                            for regex_to_substitute, regex_target in self.go_terms_replacement_dict.items():
                                sent_merger.terms = [re.sub(regex_to_substitute, regex_target, go_term_name) for
                                                     go_term_name in sent_merger.terms]
                                sent_merger.terms_ids_dict = {re.sub(regex_to_substitute, regex_target, key): value for
                                                              key, value in sent_merger.terms_ids_dict.items()}
            sentences = [GOSentence(prefix=prefix, terms=list(sent_merger.terms),
                                    term_ids_dict=sent_merger.terms_ids_dict,
                                    postfix=GOSentencesCollection.merge_postfix_phrases(sent_merger.postfix_list),
                                    text=compose_go_sentence(prefix=prefix,
                                                             go_term_names=list(sent_merger.terms),
                                                             postfix=GOSentencesCollection.merge_postfix_phrases(
                                                                 sent_merger.postfix_list),
                                                             additional_prefix=sent_merger.additional_prefix),
                                    go_aspect=go_aspect, evidence_group=", ".join(sent_merger.evidence_groups),
                                    terms_merged=True, additional_prefix=sent_merger.additional_prefix)
                         for prefix, sent_merger in merged_sentences.items() if len(sent_merger.terms) > 0]
        if desc_stats:
            desc_stats.num_terms_trim_group_priority_merge[go_aspect] = 0
            desc_stats.terms_trim_group_priority_merge[go_aspect] = []
            for sentence in sentences:
                desc_stats.num_terms_trim_group_priority_merge[sentence.go_aspect] += len(sentence.terms)
                desc_stats.terms_trim_group_priority_merge[sentence.go_aspect].extend(sentence.terms)
        return sentences

    @staticmethod
    def merge_postfix_phrases(postfix_phrases: List[str]) -> str:
        """merge postfix phrases and remove possible redundant text at the beginning at at the end of the phrases

        :param postfix_phrases: the phrases to merge
        :type postfix_phrases: List[str]
        :return: the merged postfix phrase
        :rtype: str
        """
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


def generate_go_sentences(go_annotations: List[dict], go_ontology, evidence_groups_priority_list: List[str],
                          go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]],
                          go_prepostfix_special_cases_sent_map: Dict[Tuple[str, str], Tuple[int, str, str, str]],
                          evidence_codes_groups_map: Dict[str, str], remove_parent_terms: bool = True,
                          merge_num_terms_threshold: int = 3,
                          merge_min_distance_from_root: dict = None, merge_max_distance_from_leaf: dict = None,
                          desc_stats: SingleDescStats = None,
                          go_terms_replacement_dict: Dict[str, str] = None) -> GOSentencesCollection:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[dict]
    :param go_ontology: the go ontology
    :param evidence_groups_priority_list: the list of evidence groups to consider, sorted by priority. Sentences of the
        first group (with highest priority) will be returned in first position and so on
    :type evidence_groups_priority_list: List[str]
    :param go_prepostfix_sentences_map: a map with prefix and postfix phrases, where keys are tuples of
        go_aspect, evidence_group and values are tuples prefix, postfix
    :type go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]]
    :param go_prepostfix_special_cases_sent_map: a map for special prefix and postfix cases, where keys are tuples of
        go_aspect, evidence_group and values are tuples of id, match_regex, prefix, postfix. Match_regex is a regular
        expression that defines the match for the special case
    :type go_prepostfix_special_cases_sent_map: Dict[Tuple[str, str], Tuple[int, str, str, str]]
    :param evidence_codes_groups_map: a map between evidence codes and the groups they belong to
    :type evidence_codes_groups_map: Dict[str, str]
    :param remove_parent_terms: whether to remove parent terms from the list of terms in each sentence if at least
        one children term is present
    :type remove_parent_terms: bool
    :param merge_num_terms_threshold: whether to merge terms by common ancestor to
        reduce the number of terms in the set. The trimming algorithm will be applied only if the number of terms is
        greater than the specified number and the specified threshold is greater than 0
    :type merge_num_terms_threshold: int
    :param merge_min_distance_from_root: minimum distance from root terms for the selection of common ancestors
        during merging operations. Three values must be provided in the form of a dictionary with keys 'F', 'P', and
        'C' for go aspect names and values integers indicating the threshold for each aspect
    :type merge_min_distance_from_root: dict
    :param merge_max_distance_from_leaf: maximum distance from leaf terms for the selection of common ancestors
        during merging operations. Three values must be provided in the form of a dictionary with keys 'F', 'P', and
        'C' for go aspect names and values integers indicating the threshold for each aspect
    :type merge_min_distance_from_root: dict
    :param desc_stats: an object containing the description statistics where to save the total number of annotations
        for the gene
    :type desc_stats: SingleDescStats
    :param go_terms_replacement_dict: replacement dictionary for terms
    :type go_terms_replacement_dict: Dict[str, str]
    :return: a collection of GO sentences
    :rtype: GOSentencesCollection
    """
    if len(go_annotations) > 0:
        if not merge_min_distance_from_root:
            merge_min_distance_from_root = {'F': 1, 'P': 1, 'C': 2}
        if not merge_max_distance_from_leaf:
            merge_max_distance_from_leaf = {'F': None, 'P': None, 'C': None}
        go_terms_groups = defaultdict(set)
        for annotation in go_annotations:
            if annotation["Evidence"] in evidence_codes_groups_map:
                map_key = (annotation["Aspect"], evidence_codes_groups_map[annotation["Evidence"]])
                if map_key in go_prepostfix_special_cases_sent_map:
                    for special_case in go_prepostfix_special_cases_sent_map[map_key]:
                        if re.match(re.escape(special_case[1]), annotation["GO_Name"]):
                            map_key = (annotation["Aspect"], evidence_codes_groups_map[annotation["Evidence"]] +
                                       str(special_case[0]))
                            if evidence_codes_groups_map[annotation["Evidence"]] + str(special_case[0]) not in \
                                    evidence_groups_priority_list:
                                evidence_groups_priority_list.insert(evidence_groups_priority_list.index(
                                    evidence_codes_groups_map[annotation["Evidence"]]) + 1,
                                                                     evidence_codes_groups_map[annotation["Evidence"]] +
                                                                     str(special_case[0]))
                            break
                go_terms_groups[map_key].add((annotation["GO_Name"], annotation["GO_ID"]))
        sentences = GOSentencesCollection(evidence_groups_priority_list, go_prepostfix_sentences_map,
                                          remove_parent_terms=remove_parent_terms, go_ontology=go_ontology,
                                          go_terms_replacement_dict=go_terms_replacement_dict)
        for ((go_aspect, evidence_group), go_terms) in go_terms_groups.items():
            go_term_names = [term[0] for term in go_terms]
            if desc_stats:
                desc_stats.num_terms_notrim_nogroup_priority_nomerge[go_aspect] += len(go_term_names)
                desc_stats.terms_notrim_nogroup_priority_nomerge[go_aspect].extend(go_term_names)
            term_ids_dict = {term_name: term_id for term_name, term_id in go_terms}
            replace = False
            add_others = False
            if remove_parent_terms:
                term_ids_no_parents = get_term_ids_without_parents_from_terms_names(go_terms_names=go_term_names,
                                                                                    term_ids_dict=term_ids_dict,
                                                                                    ontology=go_ontology)
                if len(go_term_names) > len(term_ids_no_parents):
                    logging.debug("Removed " + str(len(go_term_names) - len(term_ids_no_parents)) +
                                  " parents from terms")
                    go_term_names = [go_ontology.query_term(term_id).name for term_id in term_ids_no_parents]
                    term_ids_dict = {go_ontology.query_term(term_id).name: go_ontology.query_term(term_id).id for
                                     term_id in term_ids_no_parents}
                    replace = True
            if merge_num_terms_threshold > 0:
                merged_ids_covered_subsets = get_merged_term_ids_by_common_ancestor_from_term_names(
                    go_terms_names=go_term_names, term_ids_dict=term_ids_dict, ontology=go_ontology,
                    min_distance_from_root=merge_min_distance_from_root[go_aspect],
                    max_distance_from_leaf=merge_max_distance_from_leaf[go_aspect],
                    min_number_of_terms=merge_num_terms_threshold)
                if len(merged_ids_covered_subsets.keys()) <= merge_num_terms_threshold:
                    merged_ids = list(merged_ids_covered_subsets.keys())
                else:
                    merged_ids = find_set_covering([(k, v) for k, v in merged_ids_covered_subsets.items()],
                                                   max_num_subsets=merge_num_terms_threshold)
                    add_others = True
                if 0 < len(merged_ids) < len(go_term_names):
                    logging.debug("Reduced number of terms by merging from " + str(len(go_term_names)) + " to " +
                                  str(len(merged_ids)))
                    go_term_names = [go_ontology.query_term(term_id).name for term_id in merged_ids]
                    term_ids_dict = {go_ontology.query_term(term).name: go_ontology.query_term(term).id for term in
                                     merged_ids}
                    replace = True
            if go_terms_replacement_dict and replace:
                for regex_to_substitute, regex_target in go_terms_replacement_dict.items():
                    go_term_names = [re.sub(regex_to_substitute, regex_target, go_ontology.query_term(term_id).name) for
                                     term_id in merged_ids]
                    term_ids_dict = {re.sub(regex_to_substitute, regex_target, go_ontology.query_term(term).name):
                                     go_ontology.query_term(term).id for term in merged_ids}
            sentences.set_sentence(_get_single_go_sentence(go_term_names=go_term_names,
                                                           go_term_ids_dict=term_ids_dict,
                                                           go_aspect=go_aspect,
                                                           evidence_group=evidence_group,
                                                           go_prepostfix_sentences_map=go_prepostfix_sentences_map,
                                                           terms_merged=True if 0 < merge_num_terms_threshold < len(
                                                               go_term_names) else False, add_others=add_others))
        return sentences


def compose_go_sentence(prefix: str, additional_prefix: str, go_term_names: List[str], postfix: str) -> str:
    """compose the text of a sentence given its prefix, terms, and postfix

    :param prefix: the prefix of the sentence
    :type prefix: str
    :param go_term_names: a list of go terms
    :type go_term_names: List[str]
    :param postfix: the postfix of the sentence
    :type postfix: str
    :return: the text of the go sentence
    :rtype: str"""
    prefix = prefix + additional_prefix + " "
    if postfix != "":
        postfix = " " + postfix
    if len(go_term_names) > 2:
        return prefix + ", ".join(go_term_names[0:-1]) + ", and " + go_term_names[len(go_term_names) - 1] + postfix
    elif len(go_term_names) > 1:
        return prefix + " and ".join(go_term_names) + postfix
    else:
        return prefix + go_term_names[0] + postfix


def _get_single_go_sentence(go_term_names: List[str], go_term_ids_dict: Dict[str, str], go_aspect: str,
                            evidence_group: str,
                            go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]],
                            terms_merged: bool = False, add_others: bool = False) -> Union[GOSentence, None]:
    """build a go sentence

    :param go_term_names: list of go term names to be combined in the sentence
    :type go_term_names: List[str]
    :param go_term_ids_dict: map between term names and their ids
    :type go_term_ids_dict: Dict[str, str]
    :param go_aspect: go aspect
    :type go_aspect: str
    :param evidence_group: evidence group
    :type evidence_group: str
    :param go_prepostfix_sentences_map: map for prefix and postfix phrases
    :type go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]]
    :param terms_merged: whether the terms set has been merged to reduce its size
    :type terms_merged: bool
    :param add_others: whether to say that there are other terms which have been omitted from the sentence
    :type add_others: bool
    :return: the combined go sentence
    :rtype: Union[GOSentence, None]
    """
    if len(go_term_names) > 0:
        prefix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][0]
        additional_prefix = ""
        others_word = "entities"
        if go_aspect == "F":
            others_word = "functions"
        elif go_aspect == "P":
            others_word = "processes"
        elif go_aspect == "C":
            others_word = "cellular components"
        if add_others:
            additional_prefix += " several " + others_word + ", including"
        elif go_aspect == "C":
            additional_prefix += " the"
        postfix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][1]
        return GOSentence(prefix=prefix, terms=go_term_names, postfix=postfix,
                          term_ids_dict=go_term_ids_dict, text=compose_go_sentence(prefix=prefix,
                                                                                   go_term_names=go_term_names,
                                                                                   postfix=postfix,
                                                                                   additional_prefix=additional_prefix),
                          go_aspect=go_aspect, evidence_group=evidence_group, terms_merged=terms_merged,
                          additional_prefix=additional_prefix)
    else:
        return None
