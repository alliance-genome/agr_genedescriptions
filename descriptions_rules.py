from collections import namedtuple, defaultdict
from enum import Enum
from typing import List, Dict, Tuple

GOSentence = namedtuple('GOSentence', ['prefix', 'terms', 'postfix', 'text', 'go_aspect', 'evidence_group'])


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
        self.sentences_map[(sentence.go_aspect, sentence.evidence_group)] = sentence

    def get_sentences(self, go_aspect: str, keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False) -> List[GOSentence]:
        """get all sentences containing the specified aspect

        :param go_aspect: a GO aspect
        :type go_aspect: str
        :param keep_only_best_group: whether to get only the evidence group with highest priority and discard
            the other evidence groups
        :type keep_only_best_group: bool
        :param merge_groups_with_same_prefix: whether to merge the phrases for evidence groups with the same prefix
        :type merge_groups_with_same_prefix: bool
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[GOSentence]
        """
        sentences = []
        merged_sentences = defaultdict(dict)
        for eg in self.evidence_groups_list:
            if (go_aspect, eg) in self.sentences_map:
                if merge_groups_with_same_prefix:
                    if "postfix" not in merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]:
                        merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["postfix"] = []
                    merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["postfix"].append(
                        self.go_prepostfix_sentences_map[(go_aspect, eg)][1])
                    if "terms" not in merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]:
                        merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["terms"] = set()
                    merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["terms"].update(
                        self.sentences_map[(go_aspect, eg)].terms)
                    if "evidence_groups" not in merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]:
                        merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["evidence_groups"] = []
                    merged_sentences[self.go_prepostfix_sentences_map[(go_aspect, eg)][0]]["evidence_groups"].append(eg)
                else:
                    sentences.append(self.sentences_map[(go_aspect, eg)])
                if keep_only_best_group:
                    break
        if merge_groups_with_same_prefix:
            sentences = [GOSentence(prefix=prefix, terms=list(postfixes_terms["terms"]),
                                    postfix=" and ".join(postfixes_terms["postfix"]),
                                    text=compose_go_sentence(prefix=prefix, go_term_names=list(postfixes_terms["terms"]),
                                                             postfix=" and ".join(postfixes_terms["postfix"])),
                                    go_aspect=go_aspect, evidence_group=", ".join(postfixes_terms["evidence_groups"]))
                         for prefix, postfixes_terms in merged_sentences.items()]
        return sentences


def generate_go_sentences(go_annotations: List[dict], evidence_groups_priority_list: List[str],
                          go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]],
                          evidence_codes_groups_map: Dict[str, str]) -> GOSentencesCollection:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[dict]
    :param evidence_groups_priority_list: the list of evidence groups to consider, sorted by priority. Sentences of the
        first group (with highest priority) will be returned in first position and so on
    :type evidence_groups_priority_list: List[str]
    :param go_prepostfix_sentences_map: a map with for the prefix and postfix, where keys are tuples of
        go_aspect, evidence_group and values are tuples prefix, postfix
    :type go_prepostfix_sentences_map: Dict[Tuple[str, str], Tuple[str, str]]
    :param evidence_codes_groups_map: a map between evidence codes and the groups they belong to
    :type evidence_codes_groups_map: Dict[str, str]
    :return: a collection of GO sentences
    :rtype: GOSentencesCollection
    """
    if len(go_annotations) > 0:
        go_terms_groups = defaultdict(list)
        for annotation in go_annotations:
            go_terms_groups[(annotation["Aspect"], evidence_codes_groups_map[annotation["Evidence"]])].append(
                annotation["GO_Name"])
        sentences = GOSentencesCollection(evidence_groups_priority_list, go_prepostfix_sentences_map)
        for ((go_aspect, evidence_group), go_terms) in go_terms_groups.items():
            sentences.set_sentence(_get_single_go_sentence(go_term_names=go_terms, go_aspect=go_aspect,
                                                           evidence_group=evidence_group,
                                                           go_prepostfix_sentences_map=go_prepostfix_sentences_map))
        return sentences


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
    if len(go_term_names) > 1:
        return prefix + ", ".join(go_term_names[0:-1]) + ", and " + go_term_names[len(go_term_names) - 1] + postfix
    else:
        return prefix + go_term_names[0] + postfix


def _get_single_go_sentence(go_term_names: List[str], go_aspect: int, evidence_group: int,
                            go_prepostfix_sentences_map: dict):
    if len(go_term_names) > 0:
        prefix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][0]
        postfix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][1]
        return GOSentence(prefix=prefix, terms=go_term_names, postfix=postfix,
                          text=compose_go_sentence(prefix, go_term_names, postfix),
                          go_aspect=go_aspect, evidence_group=evidence_group)
