from typing import List
from enum import Enum


class GO_ASPECT(Enum):
    MOLECULAR_FUNCTION = 0
    BIOLOGICAL_PROCESS = 1
    CELLULAR_COMPONENT = 2


GO_HEADER_SENTENCES = {GO_ASPECT.MOLECULAR_FUNCTION: "encodes product that exhibits",
                       GO_ASPECT.BIOLOGICAL_PROCESS: "is involved in",
                       GO_ASPECT.CELLULAR_COMPONENT: "localizes to"}


def generate_go_sentence(functions: List[str], go_type: GO_ASPECT) -> str:
    """generate GO sentences from a list of GO terms names

    :param functions: the list of GO term names linked to a gene through a go annotation
    :type functions: List[str]
    :param go_type: type of GO function
    :type go_type: GO_ASPECT
    :return: a string containing the go sentence
    :rtype: str
    """
    if len(functions) > 0:
        if len(functions) > 1:
            return GO_HEADER_SENTENCES[go_type] + " " + ", ".join(functions[0:-1]) + ", and " + functions[len(functions) - 1]
        else:
            return GO_HEADER_SENTENCES[go_type] + " " + functions[0]
