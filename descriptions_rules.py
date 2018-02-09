from collections import namedtuple, defaultdict
from enum import Enum
from typing import List
from data_fetcher import GO_ASPECT, GOAnnotation


GOSentence = namedtuple('GOSentence', ['text', 'go_aspect', 'evidence_group'])


class GOSentencesCollection(object):
    """a group of GO sentences indexed by aspect"""
    def __init__(self, evidence_groups_list):
        self.evicdence_groups_list = evidence_groups_list
        self.sentences_map = {}

    def set_sentence(self, sentence: GOSentence) -> None:
        """add a sentence to the collection

        :param sentence: the sentence to add
        :type sentence: Sentence
        """
        self.sentences_map[(sentence.go_aspect, sentence.evidence_group)] = sentence

    def get_sentences(self, go_aspect: GO_ASPECT) -> List[GOSentence]:
        """get all sentences containing the specified aspect

        :param go_aspect: a GO aspect
        :type go_aspect: GO_ASPECT
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[Sentence]
        """
        sentences = []
        for eg in self.evicdence_groups_list:
            if (go_aspect, eg) in self.sentences_map:
                sentences.append(self.sentences_map[(go_aspect, eg)])
        return sentences


def generate_go_sentences(go_annotations: List[GOAnnotation], evidence_groups:dict, go_prepostfix_sentences_map: dict,
                          evidence_codes_groups_map: dict) -> GOSentencesCollection:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[GOAnnotation]
    :return: a collection of GO sentences
    :rtype: GOSentencesCollection
    """
    if len(go_annotations) > 0:
        go_terms_groups = defaultdict(list)
        for annotation in go_annotations:
            go_terms_groups[(annotation.aspect, evidence_codes_groups_map[annotation.evidence_code])].append(
                annotation.go_name)
        sentences = GOSentencesCollection(evidence_groups.values())
        for ((go_aspect, evidence_group), go_terms) in go_terms_groups.items():
            sentences.set_sentence(_get_single_go_sentence(go_term_names=go_terms, go_aspect=go_aspect,
                                                           evidence_group=evidence_group,
                                                           go_prepostfix_sentences_map=go_prepostfix_sentences_map))
        return sentences


def _get_single_go_sentence(go_term_names: List[str], go_aspect: int, evidence_group: int,
                            go_prepostfix_sentences_map: dict):
    if len(go_term_names) > 0:
        prefix = go_prepostfix_sentences_map[(go_aspect, evidence_group)][0] + " "
        if go_prepostfix_sentences_map[(go_aspect, evidence_group)][1] == "":
            postfix = ""
        else:
            postfix = " " + go_prepostfix_sentences_map[(go_aspect, evidence_group)][1]
        if len(go_term_names) > 1:
            return GOSentence(text=prefix + ", ".join(go_term_names[0:-1]) + ", and " +
                              go_term_names[len(go_term_names) - 1] + postfix,
                              go_aspect=go_aspect, evidence_group=evidence_group)
        else:
            return GOSentence(text=prefix + go_term_names[0] + postfix, go_aspect=go_aspect,
                              evidence_group=evidence_group)
