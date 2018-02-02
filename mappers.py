from collections import namedtuple
from enum import Enum

from namedlist import namedlist


class GO_ASPECT(Enum):
    MOLECULAR_FUNCTION = 0
    BIOLOGICAL_PROCESS = 1
    CELLULAR_COMPONENT = 2


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

Gene = namedtuple('Gene', ['id', 'name'])
GOAnnotation = namedtuple('GOAnnotation', ['go_id', 'qualifier', 'paper_reference', 'evidence_code', 'aspect',
                                           'annotation_ext', 'go_name', 'is_obsolete'])
GOOntologyEntry = namedlist('GOOntologyEntry', 'name is_obsolete')
