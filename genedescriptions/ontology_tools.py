"""set of functions to manipulate ontology graphs"""
import logging
import math
from collections import defaultdict
from typing import List, Union

from ontobio.assocmodel import AssociationSet
from ontobio.ontol import Ontology

from genedescriptions.commons import CommonAncestor

logger = logging.getLogger(__name__)


def nodes_have_same_root(node_ids: List[str], ontology: Ontology) -> Union[bool, str]:
    """
    Check whether all provided nodes are connected to the same root only

    Args:
        node_ids (List[str]): List of nodes to be checked
        ontology (Ontology): the ontology to which the provided nodes belong

    Returns:
        Union[bool, str]: the ID of the common root if all nodes are connected to the same and only root,
                          False otherwise
    """
    common_root = None
    for node_id in node_ids:
        onto_node = ontology.node(node_id)
        if "meta" in onto_node and "basicPropertyValues" in onto_node["meta"]:
            for basic_prop_val in onto_node["meta"]["basicPropertyValues"]:
                if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                    if common_root and common_root != basic_prop_val["val"]:
                        return False
                    common_root = basic_prop_val["val"]
    return True


def get_all_common_ancestors(node_ids: List[str], ontology: Ontology, min_distance_from_root: int = 0,
                             nodeids_blacklist: List[str] = None):
    """
    Retrieve all common ancestors for the provided list of nodes

    Args:
        node_ids (List[str]): list of starting nodes
        ontology (Ontology): the ontology to which the provided nodes belong
        min_distance_from_root (int): minimum distance from root node
        nodeids_blacklist (List[str]): node ids to be excluded from the result

    Returns:
        List[CommonAncestor]: list of common ancestors
    """
    common_root = nodes_have_same_root(node_ids=node_ids, ontology=ontology)
    if not common_root:
        raise ValueError("Cannot get common ancestors of nodes connected to different roots")
    ancestors = defaultdict(list)
    for node_id in node_ids:
        for ancestor in ontology.ancestors(node=node_id, reflexive=True):
            onto_anc = ontology.node(ancestor)
            onto_anc_root = None
            if "meta" in onto_anc and "basicPropertyValues" in onto_anc["meta"]:
                for basic_prop_val in onto_anc["meta"]["basicPropertyValues"]:
                    if basic_prop_val["pred"] == "OIO:hasOBONamespace":
                        onto_anc_root = basic_prop_val["val"]
            if onto_anc["depth"] >= min_distance_from_root and (
                not onto_anc_root or onto_anc_root == common_root) and (not nodeids_blacklist or ancestor not in
                                                                        nodeids_blacklist):
                ancestors[ancestor].append(node_id)
    return [CommonAncestor(node_id=ancestor, node_label=ontology.label(ancestor),
                           covered_starting_nodes=set(covered_nodes)) for ancestor, covered_nodes in
            ancestors.items() if len(covered_nodes) > 1 or ancestor == covered_nodes[0]]


def set_all_depths(ontology: Ontology, relations: List[str] = None, comparison_func=max):
    for root_id in ontology.get_roots():
        if "type" not in ontology.node(root_id) or ontology.node_type(root_id) == "CLASS":
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations,
                                       comparison_func=comparison_func)
    for node_id, node_content in ontology.nodes().items():
        if "depth" not in node_content:
            node_content["depth"] = 0


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


def set_ic_ontology_struct(ontology: Ontology, relations: List[str] = None):
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
    for node_id, node_content in ontology.nodes().items():
        if "IC" not in node_content:
            node_content["IC"] = 0
            node_content["num_subsumers"] = 0
            node_content["num_leaves"] = 0
            node_content["depth"] = 0


def reset_ic_annot_freq(ontology: Ontology):
    for node_prop in ontology.nodes().values():
        if "num_annots" in node_prop:
            del node_prop["num_annots"]
        if "IC" in node_prop:
            del node_prop["IC"]


def set_ic_annot_freq(ontology: Ontology, annotations: AssociationSet):
    for node_id, node_pr in ontology.nodes().items():
        node_pr["num_annots"] = len(annotations.query(terms=[node_id]))
    tot_annots = len(annotations.subjects)
    min_annots = min([node["num_annots"] for node in ontology.nodes().values() if "num_annots" in node and
                      node["num_annots"] > 0])
    if not min_annots:
        min_annots = 1
    for node_prop in ontology.nodes().values():
        node_prop["IC"] = -math.log(node_prop["num_annots"] / tot_annots) if node_prop["num_annots"] > 0 else -math.log(
            min_annots / (tot_annots + 1))


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



