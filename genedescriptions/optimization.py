import logging
from typing import List, Tuple, Union, Set

from ontobio import Ontology

from genedescriptions.commons import CommonAncestor


logger = logging.getLogger(__name__)


def find_set_covering(subsets: List[CommonAncestor], ontology: Ontology = None, value: List[float] = None,
                      max_num_subsets: int = None) -> Union[None, List[Tuple[str, Set[str]]]]:
    """greedy algorithm to solve set covering problem on subsets of trimming candidates

    Args:
        subsets (List[Tuple[str, str, Set[str]]]): list of subsets, each of which must contain a tuple with the first
        element being the ID of the subset, the second being the name, and the third the actual set of elements
        value (List[float]): list of costs of the subsets
        max_num_subsets (int): maximum number of subsets in the final list
    Returns:
        Union[None, List[str]]: the list of IDs of the subsets that maximize coverage with respect to the elements
                                in the element universe
    """
    logger.debug("starting set covering optimization")
    elem_to_process = {subset.node_id for subset in subsets}
    if value and len(value) != len(elem_to_process):
        return None
    universe = set([e for subset in subsets for e in subset.covered_starting_nodes])
    included_elmts = set()
    included_sets = []
    while len(elem_to_process) > 0 and included_elmts != universe and (not max_num_subsets or len(included_sets) <
                                                                       max_num_subsets):
        if value:
            effect_sets = sorted([(v * len(s.covered_starting_nodes - included_elmts), s.covered_starting_nodes,
                                   s.node_label, s.node_id) for s, v in zip(subsets, value) if s.node_id in
                                  elem_to_process], key=lambda x: (- x[0], x[2]))
        else:
            effect_sets = sorted([(len(s.covered_starting_nodes - included_elmts), s.covered_starting_nodes,
                                   s.node_label, s.node_id) for s in subsets if s.node_id in elem_to_process],
                                 key=lambda x: (- x[0], x[2]))
        elem_to_process.remove(effect_sets[0][3])
        if ontology:
            for elem in included_sets:
                if effect_sets[0][3] in ontology.ancestors(elem[0]):
                    included_sets.remove(elem)
        included_elmts |= effect_sets[0][1]
        included_sets.append((effect_sets[0][3], effect_sets[0][1]))
    logger.debug("finished set covering optimization")
    return included_sets
