"""set of functions to manipulate ontology graphs"""
import logging
import math
from collections import defaultdict
from typing import List, Tuple, Union, Set
from ontobio.ontol import Ontology


logger = logging.getLogger(__name__)


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


def set_all_information_content_values(ontology: Ontology, relations: List[str] = None):
    roots = ontology.get_roots(relations=relations)
    for root_id in roots:
        if "num_subsumers" not in ontology.node(root_id):
            _set_num_subsumers_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "num_leaves" not in ontology.node(root_id):
            _set_num_leaves_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "depth" not in ontology.node(root_id):
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        _set_information_content_in_subgraph(ontology=ontology, root_id=root_id,
                                             maxleaves=ontology.node(root_id)["num_leaves"], relations=relations)


def get_all_paths_to_root(node_id: str, ontology: Ontology, min_distance_from_root: int = 0,
                          relations: List[str] = None, nodeids_blacklist: List[str] = None,
                          previous_path: Union[None, List[str]] = None, root_node = None) -> Set[Tuple[str]]:
    """get all possible paths connecting a go term to its root terms

    Args:
        node_id (str): a valid GO id for the starting term
        ontology (Ontology): the go ontology
        min_distance_from_root (int): return only terms at a specified minimum distance from root terms
        relations (List[str]): the list of relations to be used
        nodeids_blacklist (List[str]): a list of node ids to exclude from the paths
        previous_path (Union[None, List[str]]): the path to get to the current node
    Returns:
        Set[Tuple[str]]: the set of paths connecting the specified term to its root terms, each of which contains a
        sequence of terms ids
    """
    if previous_path is None:
        previous_path = []
    new_path = previous_path[:]
    if not nodeids_blacklist or node_id not in nodeids_blacklist:
        new_path.append(node_id)
    parents = [parent for parent in ontology.parents(node=node_id, relations=relations) if
               ontology.node(parent)["depth"] >= min_distance_from_root]
    parents_same_root = []
    if root_node:
        for parent in parents:
            parent_node = ontology.node(parent)
            parent_root = None
            if "meta" in parent_node and "basicPropertyValues" in parent_node["meta"]:
                for basic_prop_val in parent_node["meta"]["basicPropertyValues"]:
                    if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                        parent_root = basic_prop_val["val"]
            if parent_root and parent_root == root_node:
                parents_same_root.append(parent)
        parents = parents_same_root

    if len(parents) > 0:
        # go up the tree, following a depth first visit
        paths_to_return = set()
        for parent in parents:
            for path in get_all_paths_to_root(node_id=parent, ontology=ontology, previous_path=new_path,
                                              min_distance_from_root=min_distance_from_root, relations=relations,
                                              nodeids_blacklist=nodeids_blacklist, root_node=root_node):
                paths_to_return.add(path)
        return paths_to_return
    if len(new_path) == 0:
        return {(node_id,)}
    else:
        return {tuple(new_path)}


def get_best_nodes_lca(node_ids: List[str], ontology: Ontology, min_distance_from_root: int = 3, max_num_nodes: int = 3,
                       nodeids_blacklist: List[str] = None) -> Tuple[bool, List[Tuple[str, Set[str]]]]:
    candidates = {node_id: (node_label, covered_nodes) for node_id, node_label, covered_nodes in
                  get_all_common_ancestors(node_ids=node_ids, ontology=ontology,
                                           min_distance_from_root=min_distance_from_root,
                                           nodeids_blacklist=nodeids_blacklist)}
    cands_ids_to_process = set(candidates.keys())
    selected_cands_ids = []
    node_to_cands_map = defaultdict(list)
    for cand in cands_ids_to_process:
        for node in candidates[cand][1]:
            node_to_cands_map[node].append(cand)
    while len(cands_ids_to_process) > 0:
        cand_id = cands_ids_to_process.pop()
        comparable_cands = [(cid, cval[1]) for cid, cval in candidates.items() if cid != cand_id and all(
            [child_id in cval[1] for child_id in candidates[cand_id][1]])]
        if len(comparable_cands) > 0:
            max_len = max(map(lambda x: len(x[1]), comparable_cands))
            best_cands = [candidate for candidate in comparable_cands if len(candidate[1]) == max_len]
            if len(best_cands) > 1:
                weighted_best_cands = sorted([(ontology.node(cand[0])["depth"], cand) for cand in best_cands],
                                             key=lambda x: x[0], reverse=True)
                max_weight = max(map(lambda x: x[0], weighted_best_cands))
                best_cands = [wcand[1] for wcand in weighted_best_cands if wcand[0] == max_weight]
            else:
                max_weight = ontology.node(best_cands[0][0])["depth"]
            if len(candidates[cand_id][1]) > len(best_cands[0][1]) or \
                    (len(candidates[cand_id][1]) > len(best_cands[0][1]) and
                     ontology.node(cand_id)["depth"] > max_weight):
                best_cands = [(cand_id, candidates[cand_id][1])]
            for best_cand in best_cands:
                selected_cands_ids.append(best_cand[0])
                for node_id in candidates[best_cand[0]][1]:
                    cands_ids_to_process -= set(node_to_cands_map[node_id])
        else:
            selected_cands_ids.append(cand_id)
    if len(selected_cands_ids) <= max_num_nodes:
        return False, [(node_id, candidates[node_id][1]) for node_id in selected_cands_ids]

    else:
        best_terms = find_set_covering([(node_id, ontology.label(node_id, id_if_null=True), candidates[node_id][1]) for
                                        node_id in selected_cands_ids], max_num_subsets=max_num_nodes)
        covered_terms = set([e for best_term_label, covered_terms in best_terms for e in covered_terms])
        return covered_terms != set(node_ids), best_terms


def get_best_nodes_naive(node_ids: List[str], ontology: Ontology, min_distance_from_root: int = 3,
                         max_num_nodes: int = 3,
                         nodeids_blacklist: List[str] = None) -> Tuple[bool, List[Tuple[str, Set[str]]]]:
    """remove terms with common ancestor and keep the ancestor term instead

    Args:
        node_ids (List[str]): the list of nodes to merge by common ancestor
        min_distance_from_root (int): set a minimum distance from root terms for ancestors that can group children terms
        max_num_nodes (int): maximum number of nodes to be returned by the trimming algorithm
        nodeids_blacklist (List[str]): a list of node ids to be excluded from common ancestors list
        ontology (Ontology): the ontology
    Returns:
        Set[str]: the set of merged terms, together with the set of original terms that each of them covers
    """
    logger.debug("applying trimming through naive algorithm")
    final_terms_set = {}
    ancestor_paths = defaultdict(list)
    term_paths = defaultdict(set)
    # step 1: get all path for each term and populate data structures
    for node_id in node_ids:
        node_root = None
        node_ont = ontology.node(node_id)
        if "meta" in node_ont and "basicPropertyValues" in node_ont["meta"]:
            for basic_prop_val in node_ont["meta"]["basicPropertyValues"]:
                if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                    node_root = basic_prop_val["val"]
        paths = get_all_paths_to_root(node_id=node_id, ontology=ontology,
                                      min_distance_from_root=min_distance_from_root, relations=None,
                                      nodeids_blacklist=nodeids_blacklist, root_node=node_root)
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
    if len(list(final_terms_set.keys())) <= max_num_nodes:
        return False, [(term_label, covered_terms) for term_label, covered_terms in final_terms_set.items()]

    else:
        best_terms = find_set_covering([(k, ontology.label(k, id_if_null=True), v) for k, v in final_terms_set.items()],
                                       max_num_subsets=max_num_nodes)
        covered_terms = set([e for best_term_label, covered_terms in best_terms for e in covered_terms])
        return covered_terms != set(node_ids), best_terms


def get_best_nodes_ic(node_ids: List[str], ontology: Ontology, max_number_of_terms: int = 3,
                      min_distance_from_root: int = 0, slim_terms_ic_bonus_perc: int = 0, slim_set: set = None,
                      nodeids_blacklist: List[str] = None) -> Tuple[bool, List[Tuple[str, Set[str]]]]:
    """trim the list of terms by selecting the best combination of terms from the initial list or their common
    ancestors based on information content

    Args:
        node_ids (List[str]): the list of nodes to merge by common ancestor
        max_number_of_terms (int): minimum number of terms above which the merge operation is performed
        ontology (Ontology): the ontology
        min_distance_from_root (int): consider only nodes at a minimum distance from root as potential candidate for
            trimming
        slim_terms_ic_bonus_perc (int): boost the IC value for terms that appear in the slim set by the provided
            percentage
        slim_set (set): set of terms that belong to the slim for the provided ontology
        nodeids_blacklist (List[str]): a list of node ids to be excluded from common ancestors list
    Returns:
        Set[str]: the set of trimmed terms, together with the set of original terms that each of them covers
    """
    common_ancestors = get_all_common_ancestors(node_ids=node_ids, ontology=ontology,
                                                nodeids_blacklist=nodeids_blacklist)
    if "IC" not in ontology.node(common_ancestors[0][0]):
        logger.warning("ontology terms do not have information content values set")
        set_all_information_content_values(ontology=ontology)
    values = [0 if node[0] not in node_ids and ontology.node(node[0])["depth"] < min_distance_from_root else
              ontology.node(node[0])["IC"] * (1 + slim_terms_ic_bonus_perc) if slim_set and node[0] in slim_set else
              ontology.node(node[0])["IC"] for node in common_ancestors]
    if slim_set and any([node[0] in slim_set for node in common_ancestors]):
        logger.debug("some candidates are present in the slim set")
    # remove ancestors with zero IC
    common_ancestors = [common_ancestor for common_ancestor, value in zip(common_ancestors, values) if value > 0]
    values = [value for value in values if value > 0]
    best_terms = find_set_covering(subsets=common_ancestors, max_num_subsets=max_number_of_terms,
                                   value=values, ontology=ontology)
    covered_terms = set([e for best_term_label, covered_terms in best_terms for e in covered_terms])
    return covered_terms != set(node_ids), best_terms


def get_all_common_ancestors(node_ids: List[str], ontology: Ontology, min_distance_from_root: int = 0,
                             nodeids_blacklist: List[str] = None):
    # check if all ids are connected to the same root node
    common_root = None
    for node_id in node_ids:
        onto_node = ontology.node(node_id)
        if "meta" in onto_node and "basicPropertyValues" in onto_node["meta"]:
            for basic_prop_val in onto_node["meta"]["basicPropertyValues"]:
                if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                    if common_root and common_root != basic_prop_val["val"]:
                        raise ValueError("Cannot get common ancestors of nodes connected to different roots")
                    common_root = basic_prop_val["val"]
    ancestors = defaultdict(list)
    for node_id in node_ids:
        for ancestor in ontology.ancestors(node=node_id, reflexive=True):
            onto_anc = ontology.node(ancestor)
            onto_anc_root = None
            if "meta" in onto_anc and "basicPropertyValues" in onto_anc["meta"]:
                for basic_prop_val in onto_anc["meta"]["basicPropertyValues"]:
                    if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                        onto_anc_root = basic_prop_val["val"]
            if onto_anc["depth"] >= min_distance_from_root and (not onto_anc_root or onto_anc_root == common_root) \
                and (not nodeids_blacklist or ancestor not in nodeids_blacklist):
                ancestors[ancestor].append(node_id)
    return [(ancestor, ontology.label(ancestor), set(covered_nodes)) for ancestor, covered_nodes in ancestors.items() if
            len(covered_nodes) > 1 or ancestor == covered_nodes[0]]


def get_best_nodes(terms, trimming_algorithm, max_terms, ontology, terms_already_covered,
                   ancestors_covering_multiple_children: Set = None, slim_bonus_perc: int = None,
                   min_dist_from_root: int = 0, slim_set=None, nodeids_blacklist: List[str] = None):
    merged_terms_coverset = None
    add_others = False
    if trimming_algorithm == "naive":
        add_others, merged_terms_coverset = get_best_nodes_naive(
            node_ids=list(terms), ontology=ontology, min_distance_from_root=min_dist_from_root,
            nodeids_blacklist=nodeids_blacklist, max_num_nodes=max_terms)
    elif trimming_algorithm == "ic":
        add_others, merged_terms_coverset = get_best_nodes_ic(
            node_ids=list(terms), ontology=ontology, max_number_of_terms=max_terms,
            min_distance_from_root=min_dist_from_root, slim_terms_ic_bonus_perc=slim_bonus_perc, slim_set=slim_set,
            nodeids_blacklist=nodeids_blacklist)
    elif trimming_algorithm == "lca":
        add_others, merged_terms_coverset = get_best_nodes_lca(
            node_ids=list(terms), ontology=ontology, min_distance_from_root=min_dist_from_root, max_num_nodes=max_terms,
            nodeids_blacklist=nodeids_blacklist)
    if ancestors_covering_multiple_children is not None:
        ancestors_covering_multiple_children.update({ontology.label(term_id, id_if_null=True) for
                                                     term_id, covered_nodes in merged_terms_coverset if
                                                     len(covered_nodes) > 1})
    terms_already_covered.update([e for term_id, covered_nodes in merged_terms_coverset for e in covered_nodes])
    terms = [term_id for term_id, covered_nodes in merged_terms_coverset]
    return terms, add_others


def _set_num_subsumers_in_subgraph(ontology: Ontology, root_id: str, relations: List[str] = None):
    if "num_subsumers" not in ontology.node(root_id):
        parents = ontology.parents(root_id)
        if not parents or all(["set_subsumers" in ontology.node(parent) for parent in parents]):
            subsumers = {subsumer for parent in parents for subsumer in ontology.node(parent)["set_subsumers"]} | \
                        {root_id}
            ontology.node(root_id)["num_subsumers"] = len(subsumers)
            ontology.node(root_id)["set_subsumers"] = subsumers
            for child_id in ontology.children(node=root_id):
                _set_num_subsumers_in_subgraph(ontology, child_id, relations)


def _set_num_leaves_in_subgraph(ontology: Ontology, root_id: str, relations: List[str] = None):
    if "set_leaves" in ontology.node(root_id):
        return ontology.node(root_id)["set_leaves"]
    children = ontology.children(node=root_id)
    if not children:
        leaves = {root_id}
        num_leaves = 0
    else:
        leaves = {leaf for child_id in children for leaf in
                  _set_num_leaves_in_subgraph(ontology=ontology, root_id=child_id, relations=relations)}
        num_leaves = len(leaves)
    ontology.node(root_id)["num_leaves"] = num_leaves
    ontology.node(root_id)["set_leaves"] = leaves
    return leaves


def _set_information_content_in_subgraph(ontology: Ontology, root_id: str, maxleaves: int, relations: List[str] = None):
    node = ontology.node(root_id)
    if str(root_id) == root_id and "ARTIFICIAL_NODE:" in root_id:
        node["IC"] = 0
    else:
        node["IC"] = -math.log((float(node["num_leaves"]) / node["num_subsumers"] + 1) / (maxleaves + 1))
    for child_id in ontology.children(node=root_id, relations=relations):
        _set_information_content_in_subgraph(ontology=ontology, root_id=child_id, maxleaves=maxleaves,
                                             relations=relations)


def find_set_covering(subsets: List[Tuple[str, str, Set[str]]], value: List[float] = None, max_num_subsets: int = None,
                      ontology: Ontology = None) -> Union[None, List[Tuple[str, Set[str]]]]:
    """greedy algorithm to solve set covering problem

    Args:
        subsets (List[Tuple[str, str, Set[str]]]): list of subsets, each of which must contain a tuple with the first
        element being the ID of the subset, the second being the name, and the third the actual set of elements
        value (List[float]): list of costs of the subsets
        max_num_subsets (int): maximum number of subsets in the final list
        ontology (Ontology): ontology to use to remove possible parent-child relationships in the result set
    Returns:
        Union[None, List[str]]: the list of IDs of the subsets that maximize coverage with respect to the elements in
        the universe
    """
    logger.debug("starting set covering optimization")
    elem_to_process = {subset[0] for subset in subsets}
    if value and len(value) != len(elem_to_process):
        return None
    universe = set([e for subset in subsets for e in subset[2]])
    included_elmts = set()
    included_sets = []
    while len(elem_to_process) > 0 and included_elmts != universe and \
            (not max_num_subsets or len(included_sets) < max_num_subsets):
        if value:
            effect_sets = sorted([(v * len(s[2] - included_elmts), s[2], s[1], s[0]) for s, v in
                                  zip(subsets, value) if s[0] in elem_to_process],
                                 key=lambda x: (- x[0], x[2]))
        else:
            effect_sets = sorted([(len(s[2] - included_elmts), s[2], s[1], s[0]) for s in subsets if s[0] in
                                  elem_to_process], key=lambda x: (- x[0], x[2]))
        elem_to_process.remove(effect_sets[0][3])
        if ontology:
            for elem in included_sets:
                if effect_sets[0][3] in ontology.ancestors(elem[0]):
                    included_sets.remove(elem)
        included_elmts |= effect_sets[0][1]
        included_sets.append((effect_sets[0][3], effect_sets[0][1]))
    logger.debug("finished set covering optimization")
    return included_sets



