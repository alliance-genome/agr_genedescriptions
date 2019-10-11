"""set of functions to manipulate ontology graphs"""
import logging
import math
from typing import List, Set
from ontobio.ontol import Ontology

from genedescriptions.trimming import TrimmingAlgorithmLCA, TrimmingAlgorithmIC, TrimmingAlgorithmNaive

logger = logging.getLogger(__name__)


def set_all_depths(ontology: Ontology, relations: List[str] = None, comparison_func=max):
    for root_id in ontology.get_roots():
        if "type" not in ontology.node(root_id) or ontology.node_type(root_id) == "CLASS":
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations,
                                       comparison_func=comparison_func)


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
        if "num_subsumers" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                              ontology.node_type(root_id) == "CLASS"):
            _set_num_subsumers_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "num_leaves" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                           ontology.node_type(root_id) == "CLASS"):
            _set_num_leaves_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "depth" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                      ontology.node_type(root_id) == "CLASS"):
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "type" not in ontology.node(root_id) or ontology.node_type(root_id) == "CLASS":
            _set_information_content_in_subgraph(ontology=ontology, root_id=root_id,
                                                 maxleaves=ontology.node(root_id)["num_leaves"], relations=relations)


def get_best_nodes(terms, trimming_algorithm, max_terms, ontology, terms_already_covered,
                   ancestors_covering_multiple_children: Set = None, slim_bonus_perc: int = None,
                   min_dist_from_root: int = 0, slim_set=None, nodeids_blacklist: List[str] = None):
    if trimming_algorithm == "ic":
        if "IC" not in ontology.node(terms[0]):
            logger.warning("ontology terms do not have information content values set")
            set_all_information_content_values(ontology=ontology)
        tr_algo = TrimmingAlgorithmIC(ontology=ontology, min_distance_from_root=min_dist_from_root,
                                      nodeids_blacklist=nodeids_blacklist, slim_terms_ic_bonus_perc=slim_bonus_perc,
                                      slim_set=slim_set)
    elif trimming_algorithm == "lca":
        tr_algo = TrimmingAlgorithmLCA(ontology=ontology, min_distance_from_root=min_dist_from_root,
                                       nodeids_blacklist=nodeids_blacklist)
    else:
        tr_algo = TrimmingAlgorithmNaive(ontology=ontology, min_distance_from_root=min_dist_from_root,
                                         nodeids_blacklist=nodeids_blacklist)
    add_others, merged_terms_coverset = tr_algo.trim(node_ids=list(terms), max_num_nodes=max_terms)
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



