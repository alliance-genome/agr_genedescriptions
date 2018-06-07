"""set of functions to manipulate ontology graphs"""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple, Union, Set
from ontobio.ontol import Ontology


def set_all_depths_in_subgraph(ontology: Ontology, root_id: str, relations: List[str] = None, comparison_func=max,
                               current_depth: int = 0):
    """calculate and set max_depth and min_depth (maximum and minimum distances from root terms in the ontology)
    recursively for all terms in a branch of the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the branch to process
        relations (List[str]): list of relations to consider
        comparison_func: a comparison function to calculate the depth when multiple paths exist between the node and
            the root. max calculates the length of the longest path, min the one of the shortest
        current_depth (int): the current depth in the ontology
    """
    if "depth" not in ontology.node(root_id):
        ontology.node(root_id)["depth"] = current_depth
    else:
        ontology.node(root_id)["depth"] = comparison_func(ontology.node(root_id)["depth"], current_depth)
    for child_id in ontology.children(node=root_id, relations=relations):
        set_all_depths_in_subgraph(ontology=ontology, root_id=child_id, relations=relations,
                                   comparison_func=comparison_func, current_depth=current_depth + 1)


def get_all_paths_to_root(node_id: str, ontology: Ontology, min_distance_from_root: int = 0,
                          relations: List[str] = None, previous_path: Union[None, List[str]] = None) -> Set[Tuple[str]]:
    """get all possible paths connecting a go term to its root terms

    Args:
        node_id (str): a valid GO id for the starting term
        ontology (Ontology): the go ontology
        min_distance_from_root (int): return only terms at a specified minimum distance from root terms
        relations (List[str]): the list of relations to be used
        previous_path (Union[None, List[str]]): the path to get to the current node
    Returns:
        Set[Tuple[str]]: the set of paths connecting the specified term to its root terms, each of which contains a
        sequence of terms ids
    """
    if previous_path is None:
        previous_path = []
    new_path = previous_path[:]
    new_path.append(node_id)
    parents = [parent for parent in ontology.parents(node=node_id, relations=relations) if
               ontology.node(parent)["depth"] >= min_distance_from_root]
    if len(parents) > 0:
        # go up the tree, following a depth first visit
        paths_to_return = set()
        for parent in parents:
            for path in get_all_paths_to_root(node_id=parent, ontology=ontology, previous_path=new_path,
                                              min_distance_from_root=min_distance_from_root, relations=relations):
                paths_to_return.add(path)
        return paths_to_return
    if len(new_path) == 0:
        return {(node_id,)}
    else:
        return {tuple(new_path)}


def get_merged_nodes_by_common_ancestor(node_ids: List[str], ontology: Ontology, min_distance_from_root: int = 3,
                                        min_number_of_terms: int = 3) -> Dict[str, Set[str]]:
    """remove terms with common ancestor and keep the ancestor term instead

    Args:
        node_ids (List[str]): the list of nodes to merge by common ancestor
        min_distance_from_root (int): set a minimum distance from root terms for ancestors that can group children terms
        min_number_of_terms (int): minimum number of terms above which the merge operation is performed
        ontology (Ontology): the ontology
    Returns:
        Set[str]: the set of merged terms, together with the set of original terms that each of them covers
    """
    if len(node_ids) > min_number_of_terms:
        logging.debug("applying trimming through naive algorithm")
        final_terms_set = {}
        ancestor_paths = defaultdict(list)
        term_paths = defaultdict(set)
        # step 1: get all path for each term and populate data structures
        for node_id in node_ids:
            paths = get_all_paths_to_root(node_id=node_id, ontology=ontology,
                                          min_distance_from_root=min_distance_from_root, relations=None)
            for path in paths:
                term_paths[node_id].add(path)
                ancestor_paths[path[-1]].append(path)
        # step 2: merge terms and keep common ancestors
        for node_id in sorted(node_ids):
            term_paths_copy = sorted(term_paths[node_id].copy(), key=lambda x: len(x))
            while len(term_paths_copy) > 0:
                curr_path = list(term_paths_copy.pop())
                selected_highest_ancestor = curr_path.pop()
                related_paths = ancestor_paths[selected_highest_ancestor]
                if not related_paths:
                    break
                covered_nodes_set = set([related_path[0] for related_path in related_paths])
                del ancestor_paths[selected_highest_ancestor]
                if curr_path:
                    if all(map(lambda x: x[0] == curr_path[0], related_paths)):
                        selected_highest_ancestor = curr_path[0]
                    else:
                        i = -1
                        while len(curr_path) > 1:
                            i -= 1
                            curr_highest_ancestor = curr_path.pop()
                            if not all(map(lambda x: len(x) >= - i, related_paths)):
                                break
                            if all(map(lambda x: x[i] == curr_highest_ancestor, related_paths)):
                                selected_highest_ancestor = curr_highest_ancestor
                                if selected_highest_ancestor in ancestor_paths:
                                    del ancestor_paths[selected_highest_ancestor]
                                for path in related_paths:
                                    term_paths[path[0]].discard(path)
                final_terms_set[selected_highest_ancestor] = covered_nodes_set
                for path in related_paths:
                    term_paths[path[0]].discard(path)
                if len(term_paths[node_id]) > 0:
                    term_paths_copy = term_paths[node_id].copy()
                else:
                    break
        logging.debug("trimming done")
        return final_terms_set
    else:
        return {node_id: {node_id} for node_id in node_ids}


def find_set_covering(subsets: List[Tuple[str, str, Set[str]]], costs: List[float] = None,
                      max_num_subsets: int = None) -> List[str]:
    """greedy algorithm to solve set covering problem

    Args:
        subsets (List[Tuple[str, str, Set[str]]]): list of subsets, each of which must contain a tuple with the first
        element being the ID of the subset, the second being the name, and the third the actual set of elements
        costs (List[float]): list of costs of the subsets
        max_num_subsets (int): maximum number of subsets in the final list
    Returns:
        List[str]: the list of IDs of the subsets that maximize coverage with respect to the elements in the universe
    """
    logging.debug("starting set covering optimization")
    if costs and len(costs) != len(subsets) and any(map(lambda x: x <= 0, costs)):
        return None
    universe = [e for subset in subsets for e in subset[1]]
    included_elmts = set()
    included_sets = []
    while len(included_sets) < len(subsets) and included_elmts != universe and \
            (not max_num_subsets or len(included_sets) < max_num_subsets):
        if costs:
            effect_sets = sorted([(c / len(s[2] - included_elmts), s[2], s[1], s[0]) for s, c in zip(subsets, costs)],
                                 key=lambda x: (- x[0], x[2]))
        else:
            effect_sets = sorted([(len(s[2] - included_elmts), s[2], s[1], s[0]) for s in subsets],
                                 key=lambda x: (- x[0], x[2]))
        included_elmts |= effect_sets[0][1]
        included_sets.append(effect_sets[0][3])
    logging.debug("finished set covering optimization")
    return included_sets



