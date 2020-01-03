from collections import namedtuple
from enum import Enum
from typing import Set, List, Any
from dataclasses import dataclass, field


@dataclass
class Sentence:
    prefix: str
    initial_terms_ids: List[str]
    terms_ids: List[str]
    postfix: str
    text: str
    aspect: str
    evidence_group: str
    terms_merged: bool
    additional_prefix: str
    qualifier: str
    ancestors_covering_multiple_terms: Set[str]
    trimmed: bool


Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


class DataType(Enum):
    GO = 1
    DO = 2
    EXPR = 3


class Module(Enum):
    GO_FUNCTION = 1
    GO_PROCESS = 2
    GO_COMPONENT = 3
    DO_EXPERIMENTAL = 4
    ORTHOLOGY = 5
    INFO_POOR = 6
    EXPRESSION = 7
    EXPRESSION_CLUSTER_GENE = 8
    EXPRESSION_CLUSTER_ANATOMY = 9
    EXPRESSION_CLUSTER_MOLECULE = 10
    DO_BIOMARKER = 11
    DO_ORTHOLOGY = 12
    SISTER_SP = 13
    INFO_POOR_HUMAN_FUNCTION = 14
    PROTEIN_DOMAIN = 15
    GO = 16
    EXPRESSION_CLUSTER_GENEREG = 17


def get_data_type_from_module(module):
    if module == Module.DO_ORTHOLOGY or module == Module.DO_EXPERIMENTAL or module == module.DO_BIOMARKER:
        return DataType.DO
    elif module == Module.GO:
        return DataType.GO
    elif module == Module.EXPRESSION:
        return DataType.EXPR
    return None


def get_module_from_data_type(data_type: DataType):
    if data_type == DataType.GO:
        return Module.GO
    elif data_type == DataType.DO:
        return Module.DO_EXPERIMENTAL
    elif data_type == DataType.EXPR:
        return Module.EXPRESSION
    return None


@dataclass
class CommonAncestor:
    node_id: Any
    node_label: str
    covered_starting_nodes: Set[str]


@dataclass
class TrimmingResult:
    final_terms: List[str] = field(default_factory=list)
    trimming_applied: bool = False
    partial_coverage: bool = False
    multicovering_nodes: Set[str] = field(default_factory=set)
    covered_nodes: Set[str] = field(default_factory=set)
