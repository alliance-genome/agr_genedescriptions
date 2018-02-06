from collections import namedtuple, defaultdict
from enum import Enum
from typing import List
from data_fetcher import GO_ASPECT, GOAnnotation


class EVIDENCE_GROUP(Enum):
    EXPERIMENTAL = 1
    SEQUENCE_BASED_ANALYSIS = 2
    PHYLOGENETIC_ANALYSIS = 3
    ELECTRONIC = 4
    COMPUTATIONAL_ANALYSIS = 5


EVIDENCE_GROUPS_MAP = {
    "EXP": EVIDENCE_GROUP.EXPERIMENTAL, "IDA": EVIDENCE_GROUP.EXPERIMENTAL,
    "IPI": EVIDENCE_GROUP.EXPERIMENTAL, "IMP": EVIDENCE_GROUP.EXPERIMENTAL,
    "IGI": EVIDENCE_GROUP.EXPERIMENTAL, "IEP": EVIDENCE_GROUP.EXPERIMENTAL,
    "ISS": EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS, "ISO": EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS,
    "ISA": EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS, "ISM": EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS,
    "IGC": EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS, "IBA": EVIDENCE_GROUP.PHYLOGENETIC_ANALYSIS,
    "IBD": EVIDENCE_GROUP.PHYLOGENETIC_ANALYSIS, "IKR": EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS,
    "IRD": EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS, "RCA": EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS,
    "IC": EVIDENCE_GROUP.EXPERIMENTAL, "IEA": EVIDENCE_GROUP.ELECTRONIC
}

GO_PREPOSTFIX_SENTENCES_MAP = {(GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.EXPERIMENTAL): ("exhibits", ""),
                               (GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS):
                               ("is predicted to have", "based on sequence-based analysis"),
                               (GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.PHYLOGENETIC_ANALYSIS):
                               ("is predicted to have", "based on phylogenetic analysis"),
                               (GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.ELECTRONIC):
                               ("is predicted to have", "based on protein domain information"),
                               (GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS):
                               ("is predicted to have", "based on computational analysis"),
                               (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.EXPERIMENTAL):
                               ("is involved in", ""),
                               (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS):
                               ("is predicted to be involved in", "based on sequence-based analysis"),
                               (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.PHYLOGENETIC_ANALYSIS):
                               ("is predicted to be involved in", "based on phylogenetic analysis"),
                               (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.ELECTRONIC):
                               ("is predicted to be involved in", "based on protein domain information"),
                               (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS):
                               ("is predicted to be involved in", "based on computational analysis"),
                               (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.EXPERIMENTAL):
                               ("localizes to", ""),
                               (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.SEQUENCE_BASED_ANALYSIS):
                               ("is predicted to localize to", "based on sequence-based analysis"),
                               (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.PHYLOGENETIC_ANALYSIS):
                               ("is predicted to localize to", "based on phylogenetic analysis"),
                               (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.ELECTRONIC):
                               ("is predicted to localize to", "based on protein domain information"),
                               (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.COMPUTATIONAL_ANALYSIS):
                               ("is predicted to localize to", "based on computational analysis")}

GOSentence = namedtuple('GOSentence', ['text', 'go_aspect', 'evidence_group'])


class GOSentencesCollection(object):
    """a group of GO sentences indexed by aspect"""
    def __init__(self):
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
        for eg in EVIDENCE_GROUP:
            if (go_aspect, eg) in self.sentences_map:
                sentences.append(self.sentences_map[(go_aspect, eg)])
        return sentences


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
            go_terms_groups[(annotation.aspect, EVIDENCE_GROUPS_MAP[annotation.evidence_code])].append(annotation.go_name)
        sentences = GOSentencesCollection()
        for ((go_aspect, evidence_group), go_terms) in go_terms_groups.items():
            sentences.set_sentence(_get_single_go_sentence(go_term_names=go_terms, go_aspect=go_aspect,
                                                           evidence_group=evidence_group))
        return sentences


def _get_single_go_sentence(go_term_names: List[str], go_aspect: GO_ASPECT, evidence_group: EVIDENCE_GROUP):
    if len(go_term_names) > 0:
        prefix = GO_PREPOSTFIX_SENTENCES_MAP[(go_aspect, evidence_group)][0] + " "
        if GO_PREPOSTFIX_SENTENCES_MAP[(go_aspect, evidence_group)][1] == "":
            postfix = ""
        else:
            postfix = " " + GO_PREPOSTFIX_SENTENCES_MAP[(go_aspect, evidence_group)][1]
        if len(go_term_names) > 1:
            return GOSentence(text=prefix + ", ".join(go_term_names[0:-1]) + ", and " +
                              go_term_names[len(go_term_names) - 1] + postfix,
                              go_aspect=go_aspect, evidence_group=evidence_group)
        else:
            return GOSentence(text=prefix + go_term_names[0] + postfix, go_aspect=go_aspect,
                              evidence_group=evidence_group)
