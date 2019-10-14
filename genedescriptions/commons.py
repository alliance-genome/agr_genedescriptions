from collections import namedtuple
from enum import Enum
from typing import Set, List, Any
from dataclasses import dataclass, field
from namedlist import namedlist

Sentence = namedlist('Sentence', ['prefix', 'terms_ids', 'postfix', 'text', 'aspect', 'evidence_group', 'terms_merged',
                                  'additional_prefix', 'qualifier', 'ancestors_covering_multiple_terms', "trimmed"])


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
