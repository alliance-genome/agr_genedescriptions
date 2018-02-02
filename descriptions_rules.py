from collections import namedtuple, defaultdict
from typing import List
from mappers import GO_ASPECT, EVIDENCE_GROUP, GOAnnotation, EVIDENCE_GROUPS_MAP

GO_HEADER_SENTENCES = {(GO_ASPECT.MOLECULAR_FUNCTION, EVIDENCE_GROUP.EXPERIMENTS): "exhibits",
                       (GO_ASPECT.BIOLOGICAL_PROCESS, EVIDENCE_GROUP.EXPERIMENTS): "is involved in",
                       (GO_ASPECT.CELLULAR_COMPONENT, EVIDENCE_GROUP.EXPERIMENTS): "localizes to"}

Sentence = namedtuple('Sentence', ['text', 'go_aspect', 'evidence_code'])


class Sentences(object):
    def __init__(self):
        self.sentences_map = defaultdict(list)

    def add_sentence(self, sentence: Sentence):
        self.sentences_map[sentence.go_aspect].append(sentence)

    def get_sentences(self, go_aspect: GO_ASPECT):
        return self.sentences_map[go_aspect]


def generate_go_sentence(go_annotations: List[GOAnnotation]) -> Sentences:
    """generate GO sentences from a list of GO annotations

    :param go_annotations: the list of GO annotations for a given gene
    :type go_annotations: List[GOAnnotation]
    :return: the list of sentences
    :rtype: List[Sentence]
    """
    if len(go_annotations) > 0:
        go_terms_groups = defaultdict(list)
        for annotation in go_annotations:
            go_terms_groups[(annotation.aspect, annotation.evidence_code)].append(annotation.go_name)
        sentences = Sentences()
        for ((go_aspect, evidence_code), go_terms) in go_terms_groups.items():
            sentences.add_sentence(_get_single_sentence(go_term_names=go_terms, go_aspect=go_aspect,
                                                        evidence_code=evidence_code))
        return sentences


def _get_single_sentence(go_term_names: List[str], go_aspect: GO_ASPECT, evidence_code: str):
    if len(go_term_names) > 0:
        if len(go_term_names) > 1:
            ec = EVIDENCE_GROUPS_MAP[evidence_code]
            return Sentence(text=GO_HEADER_SENTENCES[(go_aspect, ec)] + " " +
                            ", ".join(go_term_names[0:-1]) + ", and " + go_term_names[len(go_term_names) - 1],
                            go_aspect=go_aspect, evidence_code=evidence_code)
        else:
            ec = EVIDENCE_GROUPS_MAP[evidence_code]
            return Sentence(text=GO_HEADER_SENTENCES[(go_aspect, ec)] + " " +
                            go_term_names[0], go_aspect=go_aspect, evidence_code=evidence_code)
