import inflect
import re

from collections import namedtuple, defaultdict
from typing import List, Dict, Tuple, Union

GOSentence = namedtuple('GOSentence', ['prefix', 'terms', 'term_ids_dict', 'postfix', 'text', 'go_aspect',
                                       'evidence_group'])


class GOSentenceMerger(object):
    def __init__(self):
        self.postfix_list = []
        self.terms = set()
        self.terms_ids_dict = {}
        self.term_postfix_dict = {}
        self.evidence_groups = []
        self.term_evgroup_dict = {}


class GOSentencesCollection(object):
    """a group of GO sentences indexed by aspect"""
    def __init__(self, evidence_groups_list, go_prepostfix_sentences_map):
        self.evidence_groups_list = evidence_groups_list
        self.go_prepostfix_sentences_map = go_prepostfix_sentences_map
        self.sentences_map = {}

    def set_sentence(self, sentence: GOSentence) -> None:
        """add a sentence to the collection

        :param sentence: the sentence to add
        :type sentence: Sentence
        """
        if sentence is not None:
            self.sentences_map[(sentence.go_aspect, sentence.evidence_group)] = sentence

    def get_sentences(self, go_aspect: str, go_ontology, keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False) -> List[GOSentence]:
        """get all sentences containing the specified aspect

        :param go_aspect: a GO aspect
        :type go_aspect: str
        :param go_ontology: the go ontology object obtained from a data fetcher
        :param keep_only_best_group: whether to get only the evidence group with highest priority and discard
            the other evidence groups
        :type keep_only_best_group: bool
        :param merge_groups_with_same_prefix: whether to merge the phrases for evidence groups with the same prefix
        :type merge_groups_with_same_prefix: bool
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[GOSentence]
        """
        sentences = []
        merged_sentences = defaultdict(GOSentenceMerger)
        for eg in self.evidence_groups_list:
            if (go_aspect, eg) in self.sentences_map:
                if merge_groups_with_same_prefix:
                    prefix = self.go_prepostfix_sentences_map[(go_aspect, eg)][0]
                    merged_sentences[prefix].postfix_list.append(self.go_prepostfix_sentences_map[(go_aspect, eg)][1])
                    merged_sentences[prefix].terms.update(self.sentences_map[(go_aspect, eg)].terms)
                    merged_sentences[prefix].terms_ids_dict.update(self.sentences_map[(go_aspect, eg)].term_ids_dict)
                    for term in self.sentences_map[(go_aspect, eg)].terms:
                        merged_sentences[prefix].term_postfix_dict[term] = self.go_prepostfix_sentences_map[
                            (go_aspect, eg)][1]
                    merged_sentences[prefix].evidence_groups.append(eg)
                    for term in self.sentences_map[(go_aspect, eg)].terms:
                        merged_sentences[prefix].term_evgroup_dict[term] = eg
                else:
                    sentences.append(self.sentences_map[(go_aspect, eg)])
                if keep_only_best_group:
                    break
        if merge_groups_with_same_prefix:
            for prefix, sent_merger in merged_sentences.items():
                new_terms_set = set(sent_merger.terms)
                new_term_ids_dict = sent_merger.terms_ids_dict.copy()
                new_term_postfix_dict = sent_merger.term_postfix_dict.copy()
                new_term_evgroup_dict = sent_merger.term_evgroup_dict.copy()
                for term_name in sent_merger.terms:
                    for parent_name in get_all_go_parent_names(sent_merger.terms_ids_dict[term_name], go_ontology):
                        new_terms_set.discard(parent_name)
                        if parent_name in new_term_ids_dict:
                            del new_term_ids_dict[parent_name]
                            del new_term_postfix_dict[parent_name]
                            del new_term_evgroup_dict[parent_name]
                sent_merger.postfix_list = [postfix for postfix in sent_merger.postfix_list if
                                            postfix in set(new_term_postfix_dict.values())]
                sent_merger.evidence_groups = [eg for eg in sent_merger.evidence_groups if
                                               eg in set(new_term_evgroup_dict.values())]
                sent_merger.terms = list(new_terms_set)
                sent_merger.terms_ids_dict = new_term_ids_dict
            sentences = [GOSentence(prefix=prefix, terms=list(sent_merger.terms),
                                    term_ids_dict=sent_merger.terms_ids_dict,
                                    postfix=GOSentencesCollection.merge_postfix_phrases(sent_merger.postfix_list),
                                    text=compose_go_sentence(prefix=prefix,
                                                             go_term_names=list(sent_merger.terms),
                                                             postfix=GOSentencesCollection.merge_postfix_phrases(
                                                                 sent_merger.postfix_list)),
                                    go_aspect=go_aspect, evidence_group=", ".join(sent_merger.evidence_groups))
                         for prefix, sent_merger in merged_sentences.items() if len(sent_merger.terms) > 0]
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
                          evidence_codes_groups_map: Dict[str, str]) -> GOSentencesCollection:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[dict]
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
    :return: a collection of GO sentences
    :rtype: GOSentencesCollection
    """
    if len(go_annotations) > 0:
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
        sentences = GOSentencesCollection(evidence_groups_priority_list, go_prepostfix_sentences_map)
        for ((go_aspect, evidence_group), go_terms) in go_terms_groups.items():
            go_term_names = set([term[0] for term in go_terms])
            term_ids_dict = {term_name: term_id for term_name, term_id in go_terms}
            for go_term in go_terms:
                for parent_name in get_all_go_parent_names(go_term[1], go_ontology):
                    go_term_names.discard(parent_name)
                    if parent_name in term_ids_dict:
                        del term_ids_dict[parent_name]
            sentences.set_sentence(_get_single_go_sentence(go_term_names=list(go_term_names),
                                                           go_term_ids_dict=term_ids_dict,
                                                           go_aspect=go_aspect,
                                                           evidence_group=evidence_group,
                                                           go_prepostfix_sentences_map=go_prepostfix_sentences_map))
        return sentences


def get_all_go_parent_names(go_id, go_ontology):
    parent_names = []
    for parent in go_ontology.query_term(go_id).parents:
        parent_names.append(parent.name)
        parent_names.extend(get_all_go_parent_names(parent.id, go_ontology))
    return parent_names


def compose_go_sentence(prefix: str, go_term_names: List[str], postfix: str) -> str:
    """compose the text of a sentence given its prefix, terms, and postfix

    :param prefix: the prefix of the sentence
    :type prefix: str
    :param go_term_names: a list of go terms
    :type go_term_names: List[str]
    :param postfix: the postfix of the sentence
    :type postfix: str
    :return: the text of the go sentence
    :rtype: str"""
    prefix = prefix + " "
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
                            go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]]) -> Union[GOSentence,
                                                                                                          None]:
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
    :return: the combined go sentence
    :rtype: Union[GOSentence, None]
    """
    if len(go_term_names) > 0:
        prefix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][0]
        postfix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][1]
        return GOSentence(prefix=prefix, terms=go_term_names, postfix=postfix,
                          term_ids_dict=go_term_ids_dict, text=compose_go_sentence(prefix, go_term_names, postfix),
                          go_aspect=go_aspect, evidence_group=evidence_group)
    else:
        return None
