import logging
import re
from collections import defaultdict
from typing import Set, List, Tuple, Dict, Union
from ontobio import Ontology

from genedescriptions.api_manager import APIManager
from genedescriptions.commons import Sentence, DataType, Module

logger = logging.getLogger("Sentence generation functions")


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
                         rename_cell: bool = False) -> Union[Sentence, None]:
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
                        ancestors_covering_multiple_terms=ancestors_with_multiple_children)
    else:
        return None


def generate_ortholog_sentence_wormbase_human(orthologs: List[List[str]], human_genes_props: Dict[str, List[str]]):
    """build orthology sentence for WormBase human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol
        human_genes_props (Dict[str, List[str]]): dictionary containing human gene properties
    Returns:
        Tuple[list, str]: the orthologs and the sentence
    """
    symbol_name_arr = []
    genes_in_families = []
    gene_symbols_wo_family = []
    if len(orthologs) > 3:
        gene_families = defaultdict(list)
        gene_symbols_wo_family = set()
        for ortholog in orthologs:
            if ortholog[0] in human_genes_props and human_genes_props[ortholog[0]] and \
                    len(human_genes_props[ortholog[0]]) == 4:
                gene_families[human_genes_props[ortholog[0]][2]].append(human_genes_props[ortholog[0]])
            else:
                gene_symbols_wo_family.add(ortholog[1])
        if len(list(gene_families.keys())) == 1:
            gene_symbols_wo_family.update(set([human_p[0] + " (" + human_p[1] + ")" for human_p in
                                               gene_families[list(gene_families.keys())[0]]]))
            gene_families = {}
        else:
            for family_symbol, human_ps in gene_families.items():
                if family_symbol == "" or len(human_ps) == 1:
                    for human_p in human_ps:
                        gene_symbols_wo_family.add(human_p[0] + " (" + human_p[1] + ")")
            gene_families = {family_name: human_ps for family_name, human_ps in gene_families.items() if
                             len(human_ps) > 1 and family_name != ""}
        gene_family_names = [human_ps[0][2] + " (" + human_ps[0][3] + ")" for human_ps in gene_families.values()]
        genes_in_families = list(set([hps[0] for gene_list in gene_families.values() for hps in gene_list]))
        gene_symbols_wo_family = list(gene_symbols_wo_family)
        if len(gene_family_names) > 3:
            gene_family_names = gene_family_names[0:3]
        if len(genes_in_families) > 3:
            genes_in_families = genes_in_families[0:3]
        if len(gene_symbols_wo_family) > 3:
            gene_symbols_wo_family = gene_symbols_wo_family[0:3]
        family_word = "family"
        if len(gene_family_names) > 1:
            family_word = "families"
        sentences_arr = []
        if len(gene_symbols_wo_family) > 0:
            sentences_arr.append("human " + concatenate_words_with_oxford_comma(gene_symbols_wo_family))
        if len(gene_family_names) > 0:
            sentences_arr.append("members of the human " + concatenate_words_with_oxford_comma(gene_family_names) +
                                 " gene " + family_word + " including " + concatenate_words_with_oxford_comma(
                genes_in_families))
        orth_sentence = "is an ortholog of " + " and ".join(sentences_arr)
    else:
        symbol_name_arr = sorted([human_genes_props[best_orth[0]][0] + " (" + human_genes_props[best_orth[0]][1] +
                                  ")" if best_orth[0] in human_genes_props and human_genes_props[best_orth[0]] else
                                  best_orth[1] for best_orth in orthologs])
        orth_sentence = "is an ortholog of human " + concatenate_words_with_oxford_comma(symbol_name_arr)

    return [gene.split(" ")[0] for gene in [*genes_in_families, *gene_symbols_wo_family, *symbol_name_arr]], \
           orth_sentence


def generate_orthology_sentence_alliance_human(orthologs: List[List[str]], excluded_orthologs: bool = False):
    """build orthology sentence for Alliance human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol, and gene_name
        excluded_orthologs (bool): whether some of the orthologs have been excluded from the final set. If true, the
            final sentence will include a prefix to specify that some orthologs have been omitted
    Returns:
        str: the orthology sentence
    """
    if len(orthologs) > 0:
        prefix = "human"
        orthologs_display = sorted(orthologs, key=lambda x: x[2])
        if excluded_orthologs or len(orthologs) > 3:
            orthologs_display = orthologs_display[0:3]
            prefix = "several human genes including"
        return "orthologous to " + prefix + " " + concatenate_words_with_oxford_comma(
            [orth[1] + " (" + orth[2] + ")" if orth[2] else orth[1] for orth in orthologs_display])
    else:
        return None


def generate_ortholog_sentence_wormbase_non_c_elegans(orthologs: List[List[str]], orthologs_sp_fullname: str,
                                                      api_manager: APIManager):
    """build orthology sentence for WormBase non-human hortologs

        Args:
            orthologs (List[str]): list of human orthologs, containing gene_id, gene_symbol
            orthologs_sp_fullname (str): full name of species from which to extract orthologs
            api_manager (APIManager): api manager to send requests to wormbase and textpresso
        Returns:
            str: the orthology sentence
        """
    orth_sentence = None
    if len(orthologs) > 0:
        fullname_arr = orthologs_sp_fullname.split(" ")
        if len(fullname_arr[0]) > 2:
            fullname_arr[0] = fullname_arr[0][0] + "."
            orthologs_sp_fullname = " ".join(fullname_arr)
        if len(orthologs) > 3:
            # sort orthologs by tpc popularity and alphabetically (if tied)
            orthologs_pop = [o_p for o_p in sorted([[ortholog, api_manager.get_textpresso_popularity(ortholog[1])] for
                                                    ortholog in orthologs], key=lambda x: (x[1], x[0][1]),
                                                   reverse=True)]
            classes_orth_pop = defaultdict(list)
            orthologs_pop_wo_class = []
            for o_p in orthologs_pop:
                gene_class = api_manager.get_gene_class(o_p[0][0])
                if gene_class:
                    classes_orth_pop[gene_class].append(o_p)
                else:
                    orthologs_pop_wo_class.append(o_p)
            if len(list(classes_orth_pop.keys())) == 1:
                orthologs_pop_wo_class.extend(classes_orth_pop[list(classes_orth_pop.keys())[0]])
                classes_orth_pop = {}
            else:
                for gene_class, orths_with_pop in classes_orth_pop.items():
                    if len(orths_with_pop) == 1:
                        orthologs_pop_wo_class.extend(orths_with_pop)
            classes_orth_pop = {gene_class: ops[0] for gene_class, ops in classes_orth_pop.items() if len(ops) > 1}
            sorted_items = [[o_p, 0] for o_p in orthologs_pop_wo_class]
            sorted_items.extend([[o_p, 1, gene_class] for gene_class, o_p in classes_orth_pop.items()])
            sorted_items.sort(key=lambda x: x[0][1], reverse=True)
            if len(sorted_items) > 3:
                sorted_items = sorted_items[0:3]
            gene_symbols_wo_class = [item[0][0][1] for item in sorted_items if item[1] == 0]
            classes_symbols = [item[2] for item in sorted_items if item[1] == 1]
            genes_symbols_in_classes = [item[0][0][1] for item in sorted_items if item[1] == 1]
            sentences_arr = []
            if len(gene_symbols_wo_class) > 0:
                sentences_arr.append(orthologs_sp_fullname + " " + concatenate_words_with_oxford_comma(
                    gene_symbols_wo_class))
            if len(classes_symbols) > 0:
                genes_symbols_in_classes_sent = concatenate_words_with_oxford_comma(genes_symbols_in_classes)
                classes_symbols_sent = concatenate_words_with_oxford_comma(classes_symbols)
                classes_word = "classes" if len(classes_symbols) > 1 else "class"
                sentences_arr.append("members of the " + orthologs_sp_fullname + " " + classes_symbols_sent +
                                     " gene " + classes_word + " including " + genes_symbols_in_classes_sent)
            orth_sentence = "is an ortholog of " + " and ".join(sentences_arr)
        else:
            # sort orthologs alphabetically
            orthologs_symbols = sorted([orth[1] for orth in orthologs])
            orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + \
                            concatenate_words_with_oxford_comma(orthologs_symbols)
    return orth_sentence


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
