import logging
import re
from collections import defaultdict
from typing import Set, List, Tuple, Dict, Union
from ontobio import Ontology

from genedescriptions.api_manager import APIManager
from genedescriptions.commons import Sentence, DataType, Module

logger = logging.getLogger(__name__)


def compose_sentence(prefix: str, additional_prefix: str, term_names: List[str], postfix: str,
                     ancestors_with_multiple_children: Set[str] = None, rename_cell: bool = False) -> str:
    """compose the text of a sentence given its prefix, terms, and postfix

    Args:
        prefix (str): the prefix of the sentence
        term_names (List[str]): a list of term names
        postfix (str): the postfix of the sentence
        additional_prefix (str): an additional prefix to be used for special cases
        ancestors_with_multiple_children (Set[str]): set containing labels of terms that cover more than one children
            term in the original set and which will appear with the label '(multiple)'
    Returns:
        str: the text of the go sentence
    """
    prefix = prefix + additional_prefix + " "
    term_names = [term_name + " (multiple)" if term_name in ancestors_with_multiple_children else term_name for
                  term_name in sorted(term_names)]

    if postfix != "":
        postfix = " " + postfix
    if rename_cell:
        if "the cell" in term_names or "the Cell" in term_names:
            if len(term_names) == 1:
                prefix = prefix[0:-3]
                term_names = ["widely"]
            else:
                if not additional_prefix:
                    prefix += "several tissues including "
                term_names = [term for term in term_names if term != "the cell" and term != "the Cell"]
    return prefix + concatenate_words_with_oxford_comma(term_names) + postfix


def _get_single_sentence(node_ids: List[str], ontology: Ontology, aspect: str, evidence_group: str, qualifier: str,
                         prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]],
                         terms_merged: bool = False, add_others: bool = False,
                         truncate_others_generic_word: str = "several",
                         truncate_others_aspect_words: Dict[str, str] = None,
                         ancestors_with_multiple_children: Set[str] = None,
                         rename_cell: bool = False, trimmed: bool = False) -> Union[Sentence, None]:
    """build a sentence object

    Args:
        node_ids (List[str]): list of ids for the terms to be combined in the sentence
        ontology (Ontology): the ontology containing the nodes
        aspect (str): aspect
        evidence_group (str): evidence group
        qualifier (str): qualifier
        prepostfix_sentences_map (Dict[Tuple[str, str, str], Tuple[str, str]]): map for prefix and postfix phrases
        terms_merged (bool): whether the terms set has been merged to reduce its size
        add_others (bool): whether to say that there are other terms which have been omitted from the sentence
        truncate_others_generic_word (str): a generic word to indicate that the set of terms reported in the sentence is
            only a subset of the original terms, e.g., 'several'
        truncate_others_aspect_words (Dict[str, str]): one word for each aspect describing the kind of terms that are
            included in the aspect
        ancestors_with_multiple_children (Set[str]): set containing labels of terms that cover more than one children
            term in the original set and which will appear with the label '(multiple)'
        rename_cell (bool): whether to rename the term 'cell'
    Returns:
        Union[Sentence,None]: the combined go sentence
    """
    if len(node_ids) > 0:
        prefix = prepostfix_sentences_map[(aspect, evidence_group, qualifier)][0]
        additional_prefix = ""
        others_word = "entities"
        if aspect in truncate_others_aspect_words:
            others_word = truncate_others_aspect_words[aspect]
        if add_others:
            additional_prefix += " " + truncate_others_generic_word + " " + others_word + ", including"
        if aspect == "C":
            additional_prefix += " the"
        postfix = prepostfix_sentences_map[(aspect, evidence_group, qualifier)][1]
        term_labels = [ontology.label(node_id, id_if_null=True) for node_id in node_ids]
        return Sentence(prefix=prefix, terms_ids=node_ids, postfix=postfix,
                        text=compose_sentence(prefix=prefix, term_names=term_labels, postfix=postfix,
                                              additional_prefix=additional_prefix,
                                              ancestors_with_multiple_children=ancestors_with_multiple_children,
                                              rename_cell=rename_cell),
                        aspect=aspect, evidence_group=evidence_group, terms_merged=terms_merged,
                        additional_prefix=additional_prefix, qualifier=qualifier,
                        ancestors_covering_multiple_terms=ancestors_with_multiple_children, trimmed=trimmed)
    else:
        return None


def rename_human_ortholog_name(ortholog_name: str):
    new_name = ortholog_name
    if " family member " in ortholog_name:
        new_name = new_name.replace(" family member ", " ")
    new_name = re.sub(r"[,]? kDa$", "", new_name)
    return new_name


def is_human_ortholog_name_valid(ortholog_name: str):
    if "human uncharacterized protein" in ortholog_name.lower():
        return False
    return True


def concatenate_words_with_oxford_comma(words: List[str]):
    """concatenate words by separating them with commas and a final oxford comma if more than two or by 'and' if two

    Args:
        words (List[str]): a list of words

    Returns:
        str: a concatenated string representing the list of words
    """
    if len(words) > 2:
        return ", ".join(words[0:-1]) + ", and " + words[-1]
    else:
        return " and ".join(words)


def get_best_human_ortholog_for_info_poor(human_orthologs, ensembl_hgnc_ids_map, evidence_codes, human_df_agr, config):
    best_orth = ""
    if len(human_orthologs) > 0:
        ev_codes_group_map = config.get_evidence_codes_groups_map(module=Module.GO)
        exp_orthologs = defaultdict(int)
        predicted_orthologs = defaultdict(int)
        for ortholog in human_orthologs:
            if ortholog[0] in ensembl_hgnc_ids_map and ensembl_hgnc_ids_map[ortholog[0]]:
                ortholog_annotations = human_df_agr.get_annotations_for_gene(
                    gene_id="RGD:" + ensembl_hgnc_ids_map[ortholog[0]], annot_type=DataType.GO,
                    priority_list=evidence_codes)
                for annotation in ortholog_annotations:
                    if annotation['aspect'] == 'F':
                        if "EXPERIMENTAL" in ev_codes_group_map[annotation["evidence"]['type']]:
                            exp_orthologs[ensembl_hgnc_ids_map[ortholog[0]]] += 1
                        else:
                            predicted_orthologs[ensembl_hgnc_ids_map[ortholog[0]]] += 1
        if len(exp_orthologs) > 0:
            best_orth = sorted([(key, value) for key, value in exp_orthologs.items()], key=lambda x: x[1],
                               reverse=True)[0][0]
        elif len(predicted_orthologs) > 0:
            best_orth = sorted([(key, value) for key, value in predicted_orthologs.items()], key=lambda x: x[1],
                               reverse=True)[0][0]
    return best_orth
