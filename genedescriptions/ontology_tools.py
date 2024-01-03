"""set of functions to manipulate ontology graphs"""
import logging
import math
import time
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
    return common_root


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
    if common_root is False:
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
            if (ancestor in node_ids or onto_anc["depth"] >= min_distance_from_root) and (
                not onto_anc_root or onto_anc_root == common_root) and (not nodeids_blacklist or ancestor not in
                                                                        nodeids_blacklist):
                ancestors[ancestor].append(node_id)
    return [CommonAncestor(node_id=ancestor, node_label=ontology.label(ancestor),
                           covered_starting_nodes=set(covered_nodes)) for ancestor, covered_nodes in
            ancestors.items() if len(covered_nodes) > 1 or ancestor == covered_nodes[0]]


def set_all_depths(ontology: Ontology, relations: List[str] = None, comparison_func=max):
    logger.info("Setting depth for all nodes")
    start_time = time.time()
    for root_id in ontology.get_roots():
        if "type" not in ontology.node(root_id) or ontology.node_type(root_id) == "CLASS":
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations,
                                       comparison_func=comparison_func)
    for node_id, node_content in ontology.nodes().items():
        if "depth" not in node_content:
            node_content["depth"] = 0
    logger.info(f"setting all depths took {time.time() - start_time} seconds")


def set_all_depths_in_subgraph(ontology: Ontology, root_id: str, relations: List[str] = None, comparison_func=max,
                               current_depth: int = 0):
    """
    Calculate and set max_depth and min_depth (maximum and minimum distances from root terms in the ontology)
    for all nodes in a branch of the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the branch to process
        relations (List[str]): list of relations to consider
        comparison_func: a comparison function to calculate the depth when multiple paths exist between the node and
            the root. max calculates the length of the longest path, min the one of the shortest
        current_depth (int): the current depth in the ontology
    """
    stack = [(root_id, current_depth)]
    while stack:
        node_id, current_depth = stack.pop()
        if "depth" not in ontology.node(node_id):
            ontology.node(node_id)["depth"] = current_depth
        else:
            ontology.node(node_id)["depth"] = comparison_func(ontology.node(node_id)["depth"], current_depth)
        children = set(ontology.children(node=node_id, relations=relations))
        children.discard(node_id)
        stack.extend([(child_id, current_depth + 1) for child_id in children])


def set_ic_ontology_struct(ontology: Ontology, relations: List[str] = None):
    logger.info("Setting information content values based on ontology structure")
    start_time = time.time()
    roots = ontology.get_roots(relations=relations)
    for root_id in roots:
        if "num_subsumers" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                              ontology.node_type(root_id) == "CLASS"):
            set_num_subsumers(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "num_leaves" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                           ontology.node_type(root_id) == "CLASS"):
            set_leaf_sets(ontology=ontology, root_id=root_id, relations=relations)
            set_num_leaves(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "depth" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                      ontology.node_type(root_id) == "CLASS"):
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id, relations=relations)
    for root_id in roots:
        if "type" not in ontology.node(root_id) or ontology.node_type(root_id) == "CLASS":
            set_information_content_in_subgraph(ontology=ontology, root_id=root_id,
                                                maxleaves=ontology.node(root_id)["num_leaves"], relations=relations)
    logger.info(f"setting information content values based on ic took {time.time() - start_time} seconds")


def set_ic_annot_freq(ontology: Ontology, annotations: AssociationSet):
    logger.info("Setting information content values based on annotation frequency")
    for node_id in ontology.nodes():
        node_prop = ontology.node(node_id)
        if "rel_annot_genes" in node_prop:
            del node_prop["rel_annot_genes"]
        if "tot_annot_genes" in node_prop:
            del node_prop["tot_annot_genes"]
        if "IC" in node_prop:
            del node_prop["IC"]
    for root_id in ontology.get_roots():
        if "depth" not in ontology.node(root_id) and ("type" not in ontology.node(root_id) or
                                                      ontology.node_type(root_id) == "CLASS"):
            set_all_depths_in_subgraph(ontology=ontology, root_id=root_id)
    node_gene_map = defaultdict(set)
    for subj, obj in annotations.associations_by_subj_obj.keys():
        node_gene_map[obj].add(subj)
    for node_id in ontology.nodes():
        node_pr = ontology.node(node_id)
        node_pr["rel_annot_genes"] = node_gene_map[node_id]
    set_tot_annots(ontology)
    for node_prop in ontology.nodes().values():
        if "tot_annot_genes" not in node_prop:
            node_prop["tot_annot_genes"] = set()
    tot_annots = len(set([gene for set_genes in node_gene_map.values() for gene in set_genes]))
    min_annots = min([len(node["tot_annot_genes"]) for node in ontology.nodes().values() if "tot_annot_genes" in node
                      and len(node["tot_annot_genes"]) > 0])
    if not min_annots:
        min_annots = 1
    for node_prop in ontology.nodes().values():
        node_prop["IC"] = -math.log(len(node_prop["tot_annot_genes"]) / tot_annots) if \
            len(node_prop["tot_annot_genes"]) > 0 else -math.log(min_annots / (tot_annots + 1))
    logger.info("Finished setting information content values")


def set_tot_annots(ontology: Ontology, relations: List[str] = None):
    """
    Calculate the total number of annotated genes in a subgraph of the ontology

    Args:
        ontology (Ontology): the ontology
        relations (List[str]): list of relations to consider

    Returns:
        Set[str]: the set of all annotated genes in the subgraph
    """
    logger.info("Setting total annotation counts")
    start_time = time.time()
    for node_id in ontology.nodes():
        if "rel_annot_genes" in ontology.node(node_id) and ontology.node(node_id)["rel_annot_genes"]:
            if "tot_annot_genes" not in ontology.node(node_id):
                ontology.node(node_id)["tot_annot_genes"] = set()
            ontology.node(node_id)["tot_annot_genes"].update(ontology.node(node_id)["rel_annot_genes"])
            for ancestor_id in ontology.ancestors(node_id):
                if "tot_annot_genes" not in ontology.node(ancestor_id):
                    ontology.node(ancestor_id)["tot_annot_genes"] = set()
                ontology.node(ancestor_id)["tot_annot_genes"].update(ontology.node(node_id)["rel_annot_genes"])
    logger.info(f"setting tot annotation counts took {time.time() - start_time} seconds")


def set_num_subsumers(ontology: Ontology, root_id: str, relations: List[str] = None):
    """
    Calculate the number of subsumers for all nodes in the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the subgraph to process
        relations (List[str]): list of relations to consider
    """
    logger.info("Setting number of subsumers")
    start_time = time.time()
    visited = set()
    stack = [(root_id, set())]
    while stack:
        node_id, subsumers = stack.pop()
        subsumers = set(subsumers)
        if node_id in visited:
            continue
        parents = set(ontology.parents(node_id))
        parents.discard(node_id)
        parents = list(parents)
        if not parents or all(["set_subsumers" in ontology.node(parent) for parent in parents]):
            subsumers |= {subsumer for parent in parents for subsumer in ontology.node(parent)["set_subsumers"]} | {node_id}
            ontology.node(node_id)["num_subsumers"] = len(subsumers)
            ontology.node(node_id)["set_subsumers"] = subsumers
            children = set(ontology.children(node=node_id))
            children.discard(node_id)
            stack.extend([(child_id, subsumers) for child_id in children])
            visited.add(node_id)
    logger.info(f"setting num subsumers took {time.time() - start_time} seconds")


def set_leaf_sets(ontology: Ontology, root_id: str, relations: List[str] = None):
    """
    Set the set of leaves for each node in the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the subgraph to process
        relations (List[str]): list of relations to consider
    """
    logger.info("Setting leaf sets")
    start_time = time.time()
    visited = set()
    stack = [root_id]
    while stack:
        node_id = stack.pop()
        if node_id in visited:
            continue
        visited.add(node_id)
        if "set_leaves" in ontology.node(node_id):
            continue
        children = set(ontology.children(node=node_id, relations=relations))
        children.discard(node_id)
        if not children:
            for ancestor in ontology.ancestors(node=node_id, relations=relations):
                if "set_leaves" not in ontology.node(ancestor):
                    ontology.node(ancestor)["set_leaves"] = set()
                ontology.node(ancestor)["set_leaves"].add(node_id)
        else:
            stack.extend([child_id for child_id in children])
    logger.info(f"setting leaf sets took {time.time() - start_time} seconds")


def set_num_leaves(ontology: Ontology, root_id: str, relations: List[str] = None):
    """
    Set the set of leaves for each node in the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the subgraph to process
        relations (List[str]): list of relations to consider
    """
    logger.info("Setting number of leaves")
    start_time = time.time()
    for node_id in ontology.nodes():
        if "set_leaves" in ontology.node(node_id):
            ontology.node(node_id)["num_leaves"] = len(ontology.node(node_id)["set_leaves"])
        else:
            ontology.node(node_id)["num_leaves"] = 0
    logger.info(f"setting num leaves took {time.time() - start_time} seconds")


def set_information_content_in_subgraph(ontology: Ontology, root_id: str, maxleaves: int, relations: List[str] = None):
    """
    Calculate the information content for a node in a subgraph of the ontology

    Args:
        ontology (Ontology): the ontology
        root_id (str): the ID of the root term of the subgraph to process
        maxleaves (int): the maximum number of leaves in the subgraph
        relations (List[str]): list of relations to consider
    """
    logger.info("Calculating IC values")
    start_time = time.time()
    visited = set()
    stack = [root_id]
    while stack:
        node_id = stack.pop()
        if node_id in visited:
            continue
        visited.add(node_id)
        node = ontology.node(node_id)
        if str(node_id) == node_id and "ARTIFICIAL_NODE:" in node_id:
            node["IC"] = 0
        else:
            if "num_leaves" in node and "num_subsumers" in node:
                node["IC"] = -math.log((float(node["num_leaves"]) / node["num_subsumers"] + 1) / (maxleaves + 1))
            else:
                logger.warning("Disconnected node: " + str(node_id))
                node["IC"] = 0
        children = set(ontology.children(node=node_id, relations=relations))
        children.discard(node_id)
        stack.extend(list(children))
    logger.info(f"calculating ic values took {time.time() - start_time} seconds")


def node_is_in_branch(ontology: Ontology, node_id: str, branch_root_ids: List[str]):
    branch_root_ids = set(branch_root_ids)
    return any([parent_id in branch_root_ids for parent_id in ontology.ancestors(node=node_id, reflexive=True)])



