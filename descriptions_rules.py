from typing import List
from enum import Enum


class GO_TYPE(Enum):
    MOLECULAR_FUNCTION = 0
    BIOLOGYCAL_PROCESS = 1
    CELLULAR_COMPONENT = 2


GO_HEADER_SENTENCES = {GO_TYPE.MOLECULAR_FUNCTION: "encodes product that exhibits",
                       GO_TYPE.BIOLOGYCAL_PROCESS: "is involved in",
                       GO_TYPE.CELLULAR_COMPONENT: "localizes to"}


def generate_go_sentence(functions: List[str], go_type: GO_TYPE) -> str:
    """generate GO sentences from a list of GO functions

    :param functions: the list of GO functions for a gene
    :type functions: List[str]
    :param go_type: type of GO function
    :type go_type: GO_TYPE
    :return: a string containing the go sentence
    :rtype: str
    """
    if len(functions) > 0:
        if len(functions) > 1:
            return GO_HEADER_SENTENCES[go_type] + " " + ", ".join(functions[0:-1]) + ", and " + functions[len(functions) - 1]
        else:
            return GO_HEADER_SENTENCES[go_type] + " " + functions[0]
