"""set of functions to manipulate ontology graphs"""
import logging
from collections import defaultdict
from typing import List, Dict, Tuple, Union, Set, Any

import goatools
import pandas as pd


def get_all_go_parent_ids(go_id: str, ontology) -> List[str]:
    """get the ids of all the ancestors of a GO term, excluding the root terms

    :param go_id: a valid GO id for the starting term
    :type go_id: str
    :param ontology: the go ontology
    :return: the list of ancestors of the term
    :rtype: List[str]
    """
    parent_ids = []
    for parent in ontology.query_term(go_id).get_parents():
        # do not return root terms
        if len(ontology.query_term(parent.id).get_parents()) > 0:
            parent_ids.append(parent.id)
            parent_ids.extend(get_all_go_parent_ids(parent.id, ontology))
    return parent_ids


def get_all_term_paths_to_root(go_id: str, ontology, min_distance_from_root: int = 0,
                               previous_path: Union[None, List[str]] = None) -> Set[Tuple[str]]:
    """get all possible paths connecting a go term to its root terms

    :param go_id: a valid GO id for the starting term
    :type go_id: str
    :param ontology: the go ontology
    :param min_distance_from_root: return only terms at a specified minimum distance from root terms
    :type min_distance_from_root: int
    :param previous_path: the path to get to the current node
    :type previous_path: Union[None, List[str]]
    :return: the set of paths connecting the specified term to its root terms, each of which contains a sequence of
        terms ids
    :rtype: Set[Tuple[str]]
    """
    if previous_path is None:
        previous_path = []
    new_path = previous_path[:]
    term_properties = ontology.query_term(go_id)
    new_path.append(term_properties.id)
    parents = [parent for parent in term_properties.get_parents() if parent.depth >= min_distance_from_root]
    if len(parents) > 0:
        # go up the tree, following a depth first visit
        paths_to_return = set()
        for parent in parents:
            for path in get_all_term_paths_to_root(go_id=parent.id, ontology=ontology, previous_path=new_path,
                                                   min_distance_from_root=min_distance_from_root):
                paths_to_return.add(path)
        return paths_to_return
    if len(new_path) == 0:
        return {(go_id,)}
    else:
        return {tuple(new_path)}


def get_term_ids_without_parents_from_terms_names(go_terms_names: List[str], term_ids_dict: Dict[str, str],
                                                  ontology) -> Set[str]:
    """remove parent terms (according to the provided go ontology) from a list of terms

    :param go_terms_names: the list of go terms from which the parents will be removed
    :type go_terms_names: List[str]
    :param term_ids_dict: a dictionary that maps term names into their GO ids
    :type term_ids_dict: Dict[str, str]
    :param ontology: the go ontology
    :return: the list of parents that have been removed from the list
    :rtype: Set[str]
    """
    go_ids_set = set([term_ids_dict[term] for term in go_terms_names])
    for go_term_name in go_terms_names:
        for parent_id in get_all_go_parent_ids(term_ids_dict[go_term_name], ontology):
            go_ids_set.discard(parent_id)
    return go_ids_set


def get_merged_term_ids_by_common_ancestor_from_term_names(go_terms_names: List[str], term_ids_dict: Dict[str, str],
                                                           ontology, min_distance_from_root: int = 3,
                                                           min_number_of_terms: int = 3) -> Dict[str, Set[str]]:
    """remove terms with common ancestor and keep the ancestor term instead

    :param go_terms_names: the list of go terms from which the parents will be removed
    :type go_terms_names: List[str]
    :param term_ids_dict: a dictionary that maps term names into their GO ids
    :type term_ids_dict: Dict[str, str]
    :param min_distance_from_root: set a minimum distance from root terms for ancestors that can group children terms
    :type min_distance_from_root: int
    :param min_number_of_terms: minimum number of terms above which the merge operation is performed
    :type min_number_of_terms: int
    :param ontology: the go ontology
    :return: the set of merged terms, together with the set of original terms that each of them covers
    :rtype: Set[str]
    """
    if len(go_terms_names) > min_number_of_terms:
        logging.debug("applying trimming through naive algorithm")
        final_terms_set = {}
        ancestor_paths = defaultdict(list)
        term_paths = defaultdict(set)
        # step 1: get all path for each term and populate data structures
        for go_term_id in [term_ids_dict[term_name] for term_name in go_terms_names]:
            paths = get_all_term_paths_to_root(go_id=go_term_id, ontology=ontology,
                                               min_distance_from_root=min_distance_from_root)
            for path in paths:
                term_paths[go_term_id].add(path)
                ancestor_paths[path[-1]].append(path)
        # step 2: merge terms and keep common ancestors
        for go_term_id in sorted([term_ids_dict[term_name] for term_name in go_terms_names]):
            term_paths_copy = sorted(term_paths[go_term_id].copy(), key=lambda x: len(x))
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
                if len(term_paths[go_term_id]) > 0:
                    term_paths_copy = term_paths[go_term_id].copy()
                else:
                    break
        logging.debug("trimming done")
        return final_terms_set
    else:
        return {term_ids_dict[term]: set(term_ids_dict[term]) for term in go_terms_names}


def find_set_covering(subsets: List[Tuple[str, Set[str]]], costs: List[float] = None, max_num_subsets: int = None) -> \
        List[str]:
    """greedy algorithm to solve set covering problem

    :param subsets: list of subsets, each of which must contain a tuple with the first element being the ID of the
        subset and the second the actual set of elements
    :type subsets: List[Tuple[str, Set[str]]]
    :param costs: list of costs of the subsets
    :type costs: List[float]
    :param max_num_subsets: maximum number of subsets in the final list
    :type max_num_subsets: int
    :return: the list of IDs of the subsets that maximize coverage with respect to the elements in the universe
    :rtype: List[str]
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
            effect_sets = sorted([(c / len(s[1] - included_elmts), s[1], s[0]) for s, c in zip(subsets, costs)],
                                 key=lambda x: x[0], reverse=True)
        else:
            effect_sets = sorted([(len(s[1] - included_elmts), s[1], s[0]) for s in subsets],
                                 key=lambda x: (x[0], x[2]), reverse=True)
        included_elmts |= effect_sets[0][1]
        included_sets.append(effect_sets[0][2])
    logging.debug("finished set covering optimization")
    return included_sets



