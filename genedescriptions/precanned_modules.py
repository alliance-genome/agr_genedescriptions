from collections import defaultdict
from typing import List, Dict

from genedescriptions.api_manager import APIManager
from genedescriptions.commons import Gene, Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription
from genedescriptions.sentence_generation_functions import concatenate_words_with_oxford_comma


def set_gene_ontology_module(dm: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription,
                             gene: Gene):
    go_sent_generator_exp = OntologySentenceGenerator(gene_id=gene.id, module=Module.GO, data_manager=dm,
                                                      config=conf_parser, limit_to_group="EXPERIMENTAL")
    go_sent_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.GO, data_manager=dm,
                                                  config=conf_parser)
    contributes_to_module_sentences = go_sent_generator.get_module_sentences(
        aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if contributes_to_module_sentences.contains_sentences():
        func_module_sentences = go_sent_generator_exp.get_module_sentences(
            aspect='F', qualifier='enables', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=func_module_sentences,
                                                                   module=Module.GO_FUNCTION)
    else:
        func_module_sentences = go_sent_generator.get_module_sentences(
            aspect='F', qualifier='enables', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=func_module_sentences, module=Module.GO_FUNCTION)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=contributes_to_module_sentences, module=Module.GO_FUNCTION)
    proc_module_sentences = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='involved_in', merge_groups_with_same_prefix=True, keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences, module=Module.GO_PROCESS)
    proc_module_sentences_a1 = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='acts_upstream_of_positive_effect', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences_a1, module=Module.GO_PROCESS)
    proc_module_sentences_a2 = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='acts_upstream_of_negative_effect', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences_a2, module=Module.GO_PROCESS)
    proc_module_sentences_a3 = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='acts_upstream_of_or_within_positive_effect', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences_a3, module=Module.GO_PROCESS)
    proc_module_sentences_a4 = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='acts_upstream_of_or_within_negative_effect', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences_a4, module=Module.GO_PROCESS)
    proc_module_sentences_a5 = go_sent_generator.get_module_sentences(
        aspect='P', qualifier='acts_upstream_of_or_within', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences_a5, module=Module.GO_PROCESS)
    colocalizes_with_module_sentences = go_sent_generator.get_module_sentences(
        aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if colocalizes_with_module_sentences.contains_sentences():
        comp_module_sentences = go_sent_generator_exp.get_module_sentences(
            aspect='C', qualifier='located_in', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
        comp_module_sentences2 = go_sent_generator_exp.get_module_sentences(
            aspect='C', qualifier='part_of', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences2, module=Module.GO_COMPONENT)
    else:
        comp_module_sentences = go_sent_generator.get_module_sentences(
            aspect='C', qualifier='located_in', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
        comp_module_sentences2 = go_sent_generator.get_module_sentences(
            aspect='C', qualifier='part_of', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences2, module=Module.GO_COMPONENT)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=colocalizes_with_module_sentences,
                                                               module=Module.GO_COMPONENT)
    gene_desc.set_or_update_initial_stats(module=Module.GO_FUNCTION, sent_generator=go_sent_generator,
                                          module_sentences=contributes_to_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.GO_FUNCTION, sent_generator=go_sent_generator,
                                          module_sentences=func_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences_a1)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences_a2)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences_a3)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences_a4)
    gene_desc.set_or_update_initial_stats(module=Module.GO_PROCESS, sent_generator=go_sent_generator,
                                          module_sentences=proc_module_sentences_a5)
    gene_desc.set_or_update_initial_stats(module=Module.GO_COMPONENT, sent_generator=go_sent_generator,
                                          module_sentences=colocalizes_with_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.GO_COMPONENT, sent_generator=go_sent_generator,
                                          module_sentences=comp_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.GO_COMPONENT, sent_generator=go_sent_generator,
                                          module_sentences=comp_module_sentences2)


def set_disease_module(df: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription, gene: Gene,
                       human: bool = False):
    do_sentence_exp_generator = OntologySentenceGenerator(gene_id=gene.id,
                                                          module=Module.DO_EXPERIMENTAL, data_manager=df,
                                                          config=conf_parser, limit_to_group="EXPERIMENTAL",
                                                          humans=human)
    disease_exp_module_sentences = do_sentence_exp_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_EXPERIMENTAL,
                                                               module_sentences=disease_exp_module_sentences)
    do_sentence_bio_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.DO_BIOMARKER,
                                                          data_manager=df, config=conf_parser,
                                                          limit_to_group="BIOMARKER", humans=human)
    disease_bio_module_sentences = do_sentence_bio_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_BIOMARKER,
                                                               module_sentences=disease_bio_module_sentences)
    do_via_orth_sentence_generator = OntologySentenceGenerator(
        gene_id=gene.id, module=Module.DO_ORTHOLOGY, data_manager=df, config=conf_parser, humans=human)
    disease_via_orth_module_sentences = do_via_orth_sentence_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_ORTHOLOGY,
                                                               module_sentences=disease_via_orth_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.DO_EXPERIMENTAL, sent_generator=do_sentence_exp_generator,
                                          module_sentences=disease_exp_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.DO_BIOMARKER, sent_generator=do_sentence_bio_generator,
                                          module_sentences=disease_bio_module_sentences)
    gene_desc.set_or_update_initial_stats(module=Module.DO_ORTHOLOGY, sent_generator=do_via_orth_sentence_generator,
                                          module_sentences=disease_via_orth_module_sentences)


def set_expression_module(df: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription, gene: Gene):
    expr_sentence_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.EXPRESSION, data_manager=df,
                                                        config=conf_parser)
    expression_module_sentences = expr_sentence_generator.get_module_sentences(
        aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=expression_module_sentences,
                                                               module=Module.EXPRESSION)
    gene_desc.set_or_update_initial_stats(module=Module.EXPRESSION, sent_generator=expr_sentence_generator,
                                          module_sentences=expression_module_sentences)


def set_alliance_human_orthology_module(orthologs: List[List[str]], gene_desc: GeneDescription,
                                        config: GenedescConfigParser, excluded_orthologs: bool = False):
    """set orthology module for Alliance human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol, and gene_name
        gene_desc (GeneDescription): the gene description object to update
        config (GenedescConfigParser): a gene descriptions configuration object
        excluded_orthologs (bool): whether some of the orthologs have been excluded from the final set. If true, the
            final sentence will include a prefix to specify that some orthologs have been omitted
    """
    if len(orthologs) > 0:
        prefix = "human"
        orthologs_display = sorted(orthologs, key=lambda x: x[2])
        if excluded_orthologs or len(orthologs) > 3:
            orthologs_display = orthologs_display[0:3]
            prefix = "several human genes including"
        sentence = "orthologous to " + prefix + " " + concatenate_words_with_oxford_comma(
            [orth[1] + " (" + orth[2] + ")" if orth[2] else orth[1] for orth in orthologs_display],
            separator=config.get_terms_delimiter())
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.ORTHOLOGY, description=sentence)


def generate_ortholog_sentence_wormbase_human(orthologs: List[List[str]], human_genes_props: Dict[str, List[str]],
                                              config: GenedescConfigParser):
    """build orthology sentence for WormBase human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol
        human_genes_props (Dict[str, List[str]]): dictionary containing human gene properties
        config (GenedescConfigParser): a gene description configuration object
    Returns:
        Tuple[list, str]: the orthologs and the sentence
    """
    prefix = "human "
    if len(orthologs) > 3:
        orthologs = orthologs[0:3]
        prefix = "several human genes including "
    symbol_name_arr = sorted([human_genes_props[best_orth[0]][0] + " (" + human_genes_props[best_orth[0]][1] +
                              ")" if best_orth[0] in human_genes_props and human_genes_props[best_orth[0]] else
                              best_orth[1] for best_orth in orthologs])
    orth_sentence = "is an ortholog of " + prefix + concatenate_words_with_oxford_comma(
        symbol_name_arr, separator=config.get_terms_delimiter())
    return [human_genes_props[best_orth[0]][0] for best_orth in orthologs if best_orth[0] in human_genes_props and
            human_genes_props[best_orth[0]]], orth_sentence


def generate_ortholog_sentence_wormbase_non_c_elegans(orthologs: List[List[str]], orthologs_sp_fullname: str,
                                                      api_manager: APIManager, config: GenedescConfigParser):
    """build orthology sentence for WormBase non-human hortologs

        Args:
            orthologs (List[str]): list of human orthologs, containing gene_id, gene_symbol
            orthologs_sp_fullname (str): full name of species from which to extract orthologs
            api_manager (APIManager): api manager to send requests to wormbase and textpresso
            config (GenedescConfigParser): a gene description configuration object
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
                    gene_symbols_wo_class, separator=config.get_terms_delimiter()))
            if len(classes_symbols) > 0:
                genes_symbols_in_classes_sent = concatenate_words_with_oxford_comma(
                    genes_symbols_in_classes, separator=config.get_terms_delimiter())
                classes_symbols_sent = concatenate_words_with_oxford_comma(classes_symbols,
                                                                           separator=config.get_terms_delimiter())
                classes_word = "classes" if len(classes_symbols) > 1 else "class"
                sentences_arr.append("members of the " + orthologs_sp_fullname + " " + classes_symbols_sent +
                                     " gene " + classes_word + " including " + genes_symbols_in_classes_sent)
            orth_sentence = "is an ortholog of " + " and ".join(sentences_arr)
        else:
            # sort orthologs alphabetically
            orthologs_symbols = sorted([orth[1] for orth in orthologs])
            orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + \
                            concatenate_words_with_oxford_comma(orthologs_symbols,
                                                                separator=config.get_terms_delimiter())
    return orth_sentence
