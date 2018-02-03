from collections import namedtuple, defaultdict
from enum import Enum
from typing import List
from data_fetcher import GO_ASPECT, GOAnnotation


class EVIDENCE_GROUP(Enum):
    EXPERIMENTS = 1


EVIDENCE_GROUPS_MAP = {
    "EXP": EVIDENCE_GROUP.EXPERIMENTS, "IDA": EVIDENCE_GROUP.EXPERIMENTS,
    "IPI": EVIDENCE_GROUP.EXPERIMENTS, "IMP": EVIDENCE_GROUP.EXPERIMENTS,
    "IGI": EVIDENCE_GROUP.EXPERIMENTS, "IEP": EVIDENCE_GROUP.EXPERIMENTS,
    "ISS": EVIDENCE_GROUP.EXPERIMENTS, "ISO": EVIDENCE_GROUP.EXPERIMENTS,
    "ISA": EVIDENCE_GROUP.EXPERIMENTS, "ISM": EVIDENCE_GROUP.EXPERIMENTS,
    "IGC": EVIDENCE_GROUP.EXPERIMENTS, "IBA": EVIDENCE_GROUP.EXPERIMENTS,
    "IBD": EVIDENCE_GROUP.EXPERIMENTS, "IKR": EVIDENCE_GROUP.EXPERIMENTS,
    "IRD": EVIDENCE_GROUP.EXPERIMENTS, "RCA": EVIDENCE_GROUP.EXPERIMENTS,
    "IC": EVIDENCE_GROUP.EXPERIMENTS
}

GO_HEADER_SENTENCES_MAP = {(GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.EXPERIMENTS): "exhibits",
                           (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.EXPERIMENTS): "is involved in",
                           (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.EXPERIMENTS): "localizes to"}

GOSentence = namedtuple('GOSentence', ['text', 'go_aspect', 'evidence_code'])


class GOSentencesCollection(object):
    """a group of GO sentences indexed by aspect"""
    def __init__(self):
        self.sentences_map = defaultdict(list)

    def add_sentence(self, sentence: GOSentence) -> None:
        """add a sentence to the collection

        :param sentence: the sentence to add
        :type sentence: Sentence
        """
        self.sentences_map[sentence.go_aspect].append(sentence)

    def get_sentences(self, go_aspect: GO_ASPECT) -> List[GOSentence]:
        """get all sentences containing the specified aspect

        :param go_aspect: a GO aspect
        :type go_aspect: GO_ASPECT
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[Sentence]
        """
        return self.sentences_map[go_aspect]


def generate_go_sentences(go_annotations: List[GOAnnotation]) -> GOSentencesCollection:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[GOAnnotation]
    :return: a collection of GO sentences
    :rtype: GOSentencesCollection
    """
    if len(go_annotations) > 0:
        go_terms_groups = defaultdict(list)
        for annotation in go_annotations:
            go_terms_groups[(annotation.aspect, annotation.evidence_code)].append(annotation.go_name)
        sentences = GOSentencesCollection()
        for ((go_aspect, evidence_code), go_terms) in go_terms_groups.items():
            sentences.add_sentence(_get_single_go_sentence(go_term_names=go_terms, go_aspect=go_aspect,
                                                           evidence_code=evidence_code))
        return sentences


def _get_single_go_sentence(go_term_names: List[str], go_aspect: GO_ASPECT, evidence_code: str):
    if len(go_term_names) > 0:
        if len(go_term_names) > 1:
            ec = EVIDENCE_GROUPS_MAP[evidence_code]
            return GOSentence(text=GO_HEADER_SENTENCES_MAP[(go_aspect, ec)] + " " + ", ".join(go_term_names[0:-1]) +
                                   ", and " + go_term_names[len(go_term_names) - 1],
                              go_aspect=go_aspect, evidence_code=evidence_code)
        else:
            ec = EVIDENCE_GROUPS_MAP[evidence_code]
            return GOSentence(text=GO_HEADER_SENTENCES_MAP[(go_aspect, ec)] + " " +
                                   go_term_names[0], go_aspect=go_aspect, evidence_code=evidence_code)
