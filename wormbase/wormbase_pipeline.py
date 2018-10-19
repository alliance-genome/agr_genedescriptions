#!/usr/bin/env python3

import argparse
import datetime

from genedescriptions.data_manager import WBDataManager, DataManager
from genedescriptions.descriptions_generator import *
from genedescriptions.descriptions_generator import generate_ortholog_sentence_wormbase_human, \
    generate_ortholog_sentence_wormbase_non_c_elegans
from genedescriptions.descriptions_writer import DescriptionsWriter


def load_data(species, organism, conf_parser):
    sister_df = None
    df_agr = None
    sister_sp_fullname = ""
    if "main_sister_species" in species[organism] and "full_name" in \
            species[species[organism]["main_sister_species"]]:
        sister_sp_fullname = species[species[organism]["main_sister_species"]]["full_name"]
    orthologs_sp_fullname = ""
    if "ortholog" in species[organism] and all(["full_name" in species[ortholog_sp] for ortholog_sp in
                                                species[organism]["ortholog"]]):
        orthologs_sp_fullname = [species[ortholog_sp]["full_name"] for ortholog_sp in species[organism]["ortholog"]]
    df = WBDataManager(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                       release_version=conf_parser.get_release("wb_data_fetcher"),
                       species=organism, project_id=species[organism]["project_id"],
                       cache_location=conf_parser.get_cache_location(), do_relations=None,
                       go_relations=["subClassOf", "BFO:0000050"], sister_sp_fullname=sister_sp_fullname)
    if organism == "c_elegans":
        df_agr = DataManager(go_relations=["subClassOf", "BFO:0000050"], do_relations=None)
        df_agr.load_ontology_from_file(ontology_type=DataType.GO,
                                       ontology_url=conf_parser.get_wb_human_orthologs_go_ontology(),
                                       ontology_cache_path=os.path.join(
                                           conf_parser.get_cache_location(), "wormbase_agr_human",
                                           "go_ontology.obo"),
                                       terms_replacement_regex=conf_parser.get_go_rename_terms())
        df_agr.load_associations_from_file(associations_type=DataType.GO,
                                           associations_url=conf_parser.get_wb_human_orthologs_go_associations(),
                                           associations_cache_path=os.path.join(
                                               conf_parser.get_cache_location(), "wormbase_agr_human",
                                               "go_assoc.daf.gz"),
                                           exclusion_list=conf_parser.get_go_terms_exclusion_list())
    if "main_sister_species" in species[organism] and species[organism]["main_sister_species"]:
        sister_df = WBDataManager(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                                  release_version=conf_parser.get_release("wb_data_fetcher"),
                                  species=species[organism]["main_sister_species"],
                                  project_id=species[species[organism]["main_sister_species"]]["project_id"],
                                  cache_location=conf_parser.get_cache_location(), do_relations=None,
                                  go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading all data for sister species")
        sister_df.load_all_data_from_file(go_terms_replacement_regex=conf_parser.get_go_rename_terms(),
                                          go_terms_exclusion_list=conf_parser.get_go_terms_exclusion_list(),
                                          do_terms_replacement_regex=None,
                                          do_terms_exclusion_list=conf_parser.get_do_terms_exclusion_list())
    logger.info("Loading all data for main species")
    df.load_all_data_from_file(go_terms_replacement_regex=conf_parser.get_go_rename_terms(),
                               go_terms_exclusion_list=conf_parser.get_go_terms_exclusion_list(),
                               do_terms_replacement_regex=None,
                               do_terms_exclusion_list=conf_parser.get_do_terms_exclusion_list())
    if organism == "c_elegans":
        df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url=df.expression_ontology_url,
                                   ontology_cache_path=df.expression_ontology_cache_path,
                                   terms_replacement_regex=conf_parser.get_expression_rename_terms())
        df.load_associations_from_file(associations_type=DataType.EXPR,
                                       associations_url=df.expression_associations_url,
                                       associations_cache_path=df.expression_associations_cache_path,
                                       exclusion_list=conf_parser.get_expression_terms_exclusion_list())
        df.load_expression_enriched_extra_info()
    elif organism == "b_malayi":
        df.load_bma_expression_data()
    elif organism == "p_pacificus":
        df.load_ppa_expression_data()
    return df, sister_df, df_agr, orthologs_sp_fullname, sister_sp_fullname


def add_orthology_sentence(df, gene, orthologs_sp_fullname, human_genes_props, gene_desc, conf_parser, joined_sent):
    best_orthologs, selected_orth_name = df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=orthologs_sp_fullname)
    selected_orthologs = []
    if best_orthologs:
        gene_desc.stats.set_best_orthologs = [orth[0] for orth in best_orthologs]
        if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
            sel_orthologs, orth_sent = generate_ortholog_sentence_wormbase_human(best_orthologs, human_genes_props)
            selected_orthologs = [orth for orth in best_orthologs if orth[1] in sel_orthologs]
        else:
            orth_sent = generate_ortholog_sentence_wormbase_non_c_elegans(best_orthologs, selected_orth_name,
                                                                          conf_parser.get_textpresso_api_token())
        if orth_sent:
            joined_sent.append(orth_sent)
            gene_desc.orthology_description = orth_sent
    return selected_orthologs


def set_initial_go_stats(go_annotations, go_sent_generator, go_sent_generator_exp, gene_desc, conf_parser):
    gene_desc.stats.total_number_go_annotations = len(go_annotations)
    gene_desc.stats.set_initial_experimental_go_ids_f = list(set().union(
        [elem for key, sets in go_sent_generator_exp.terms_groups[('F', '')].items() for elem in sets if ('F', key, '')
         in conf_parser.get_go_prepostfix_sentences_map()], [elem for key, sets in go_sent_generator_exp.terms_groups[
            ('F', 'contributes_to')].items() for elem in sets if ('F', key, 'contributes_to') in
                                                             conf_parser.get_go_prepostfix_sentences_map()]))
    gene_desc.stats.set_initial_go_ids_f = list(set().union(
        [elem for key, sets in go_sent_generator.terms_groups[('F', '')].items() for elem in sets if ('F', key, '') in
         conf_parser.get_go_prepostfix_sentences_map()], [elem for key, sets in go_sent_generator.terms_groups[
            ('F', 'contributes_to')].items() for elem in sets if ('F', key, 'contributes_to') in
                                                          conf_parser.get_go_prepostfix_sentences_map()]))
    gene_desc.stats.set_initial_go_ids_p = [elem for key, sets in go_sent_generator.terms_groups[('P', '')].items() for
                                            elem in sets if ('P', key, '') in
                                            conf_parser.get_go_prepostfix_sentences_map()]
    gene_desc.stats.set_initial_experimental_go_ids_p = [elem for key, sets in
                                                         go_sent_generator_exp.terms_groups[('P', '')].items() for
                                                         elem in sets if ('P', key, '') in
                                                         conf_parser.get_go_prepostfix_sentences_map()]
    gene_desc.stats.set_initial_go_ids_c = list(set().union(
        [elem for key, sets in go_sent_generator.terms_groups[('C', '')].items() for elem in sets if ('C', key, '') in
         conf_parser.get_go_prepostfix_sentences_map()],
        [elem for key, sets in go_sent_generator.terms_groups[('C', 'colocalizes_with')].items() for elem in sets if
         ('C', key, 'colocalizes_with') in conf_parser.get_go_prepostfix_sentences_map()]))
    gene_desc.stats.set_initial_experimental_go_ids_c = list(set().union(
        [elem for key, sets in go_sent_generator_exp.terms_groups[('C', '')].items() for elem in sets if ('C', key, '')
         in conf_parser.get_go_prepostfix_sentences_map()],
        [elem for key, sets in go_sent_generator_exp.terms_groups[('C', 'colocalizes_with')].items() for elem in sets if
         ('C', key, 'colocalizes_with') in conf_parser.get_go_prepostfix_sentences_map()]))


def set_go_function_sentence_and_set_stats(go_sent_generator, go_sent_generator_exp, go_sent_common_props, gene_desc,
                                           joined_sent):
    contributes_to_raw_func_sent = go_sent_generator.get_sentences(
        aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True, keep_only_best_group=True,
        **go_sent_common_props)
    if contributes_to_raw_func_sent:
        raw_func_sent = go_sent_generator_exp.get_sentences(aspect='F', merge_groups_with_same_prefix=True,
                                                            keep_only_best_group=True, **go_sent_common_props)
    else:
        raw_func_sent = go_sent_generator.get_sentences(aspect='F', merge_groups_with_same_prefix=True,
                                                        keep_only_best_group=True, **go_sent_common_props)
    func_sent = " and ".join([sentence.text for sentence in raw_func_sent])
    if func_sent:
        joined_sent.append(func_sent)
        gene_desc.go_function_description = func_sent
        gene_desc.go_description = func_sent
    gene_desc.stats.set_final_go_ids_f = list(set().union([term_id for sentence in raw_func_sent for
                                                           term_id in sentence.terms_ids],
                                                          [term_id for sentence in contributes_to_raw_func_sent for
                                                           term_id in sentence.terms_ids]))
    gene_desc.stats.set_final_experimental_go_ids_f = list(set().union(
        [term_id for sentence in raw_func_sent for term_id in sentence.terms_ids if
         sentence.evidence_group.startswith("EXPERIMENTAL")], [term_id for sentence in contributes_to_raw_func_sent for
                                                               term_id in sentence.terms_ids if
                                                               sentence.evidence_group.startswith("EXPERIMENTAL")]))
    contributes_to_func_sent = " and ".join([sentence.text for sentence in contributes_to_raw_func_sent])
    if contributes_to_func_sent:
        joined_sent.append(contributes_to_func_sent)
        if not gene_desc.go_function_description:
            gene_desc.go_function_description = contributes_to_func_sent
        else:
            gene_desc.go_function_description += "; " + contributes_to_func_sent
        if not gene_desc.go_description:
            gene_desc.go_description = contributes_to_func_sent
        else:
            gene_desc.go_description += "; " + contributes_to_func_sent


def set_go_process_sentence_and_set_stats(go_sent_generator, go_sent_common_props, gene_desc, joined_sent):
    raw_proc_sent = go_sent_generator.get_sentences(aspect='P', merge_groups_with_same_prefix=True,
                                                    keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.set_final_go_ids_p = [term_id for sentence in raw_proc_sent for term_id in sentence.terms_ids]
    gene_desc.stats.set_final_experimental_go_ids_p = [term_id for sentence in raw_proc_sent for term_id in
                                                       sentence.terms_ids if
                                                       sentence.evidence_group.startswith("EXPERIMENTAL")]
    proc_sent = " and ".join([sentence.text for sentence in raw_proc_sent])
    if proc_sent:
        joined_sent.append(proc_sent)
        gene_desc.go_process_description = proc_sent
        if not gene_desc.go_description:
            gene_desc.go_description = proc_sent
        else:
            gene_desc.go_description += "; " + proc_sent


def set_go_component_sentence_and_set_stats(go_sent_generator, go_sent_generator_exp, go_sent_common_props, gene_desc,
                                            joined_sent):
    colocalizes_with_raw_comp_sent = go_sent_generator.get_sentences(
        aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True, **go_sent_common_props)
    if colocalizes_with_raw_comp_sent:
        raw_comp_sent = go_sent_generator_exp.get_sentences(aspect='C', merge_groups_with_same_prefix=True,
                                                            keep_only_best_group=True, **go_sent_common_props)
    else:
        raw_comp_sent = go_sent_generator.get_sentences(aspect='C', merge_groups_with_same_prefix=True,
                                                        keep_only_best_group=True, **go_sent_common_props)
    comp_sent = " and ".join([sentence.text for sentence in raw_comp_sent])
    if comp_sent:
        joined_sent.append(comp_sent)
        gene_desc.go_component_description = comp_sent
        if not gene_desc.go_description:
            gene_desc.go_description = comp_sent
        else:
            gene_desc.go_description += "; " + comp_sent

    gene_desc.stats.set_final_go_ids_c = list(set().union([term_id for sentence in raw_comp_sent for
                                                           term_id in sentence.terms_ids],
                                                          [term_id for sentence in colocalizes_with_raw_comp_sent for
                                                           term_id in sentence.terms_ids]))
    gene_desc.stats.set_final_experimental_go_ids_c = list(set().union([term_id for sentence in raw_comp_sent for
                                                                        term_id in sentence.terms_ids if
                                                                        sentence.evidence_group.startswith(
                                                                            "EXPERIMENTAL")],
                                                                       [term_id for sentence in
                                                                        colocalizes_with_raw_comp_sent for
                                                                        term_id in sentence.terms_ids if
                                                                        sentence.evidence_group.startswith(
                                                                            "EXPERIMENTAL")]))
    colocalizes_with_comp_sent = " and ".join([sentence.text for sentence in colocalizes_with_raw_comp_sent])
    if colocalizes_with_comp_sent:
        joined_sent.append(colocalizes_with_comp_sent)
        if not gene_desc.go_component_description:
            gene_desc.go_component_description = colocalizes_with_comp_sent
        else:
            gene_desc.go_component_description += "; " + colocalizes_with_comp_sent
        if not gene_desc.go_description:
            gene_desc.go_description = colocalizes_with_comp_sent
        else:
            gene_desc.go_description += "; " + colocalizes_with_comp_sent


def set_expression_sentence(expr_sentence_generator, expr_sent_common_props, gene_desc, gene, joined_sent, df):
    raw_expression_sent = expr_sentence_generator.get_sentences(
        aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True, keep_only_best_group=False,
        **expr_sent_common_props)
    expression_sent = "; ".join([sentence.text for sentence in raw_expression_sent])
    if expression_sent:
        gene_desc.tissue_expression_description = expression_sent
        joined_sent.append(expression_sent)
        gene_desc.stats.set_final_expression_ids = [term_id for sentence in raw_expression_sent for term_id in
                                                    sentence.terms_ids]
    if len(joined_sent) == 0:
        raw_expression_sent_enriched = expr_sentence_generator.get_sentences(
            aspect='A', qualifier="Enriched", merge_groups_with_same_prefix=True, keep_only_best_group=False,
            **expr_sent_common_props)
        expression_sent_enriched = ""
        postfix = ""
        if df.expression_ontology is not None:
            expression_sent_enriched = "; ".join([sentence.text for sentence in raw_expression_sent_enriched])
            postfix = " " + concatenate_words_with_oxford_comma(
                df.expression_enriched_extra_data[gene.id[3:]]) + " studies"
            if expression_sent_enriched:
                gene_desc.gene_expression_cluster_description = expression_sent_enriched + postfix
        elif df.expression_enriched_bma_data[gene.id[3:]] and len(df.expression_enriched_bma_data[gene.id[3:]][3]) > 0:
            expression_sent_enriched = "is enriched in " + concatenate_words_with_oxford_comma(
                df.expression_enriched_bma_data[gene.id[3:]][2])
            postfix = " based on " + concatenate_words_with_oxford_comma(
                df.expression_enriched_bma_data[gene.id[3:]][3]) + " studies"
            if expression_sent_enriched:
                gene_desc.anatomy_expression_cluster_description = expression_sent_enriched + postfix
        elif df.expression_enriched_ppa_data[gene.id[3:]] and len(df.expression_enriched_ppa_data[gene.id[3:]][3]) > 0:
            expression_sent_enriched = "is enriched in " + concatenate_words_with_oxford_comma(
                df.expression_enriched_ppa_data[gene.id[3:]][2])
            postfix = " based on " + concatenate_words_with_oxford_comma(
                df.expression_enriched_ppa_data[gene.id[3:]][3]) + " studies"
            if expression_sent_enriched:
                gene_desc.anatomy_expression_cluster_description = expression_sent_enriched + postfix
        if expression_sent_enriched:
            joined_sent.append(expression_sent_enriched + postfix)
        if df.expression_ontology is None and df.expression_affected_bma_data[gene.id[3:]] and \
                len(df.expression_affected_bma_data[gene.id[3:]][3]) > 0:
            expression_sent_affected = "is affected by " + concatenate_words_with_oxford_comma(
                df.expression_affected_bma_data[gene.id[3:]][2]) + " based on " + \
                                       concatenate_words_with_oxford_comma(
                                           df.expression_affected_bma_data[gene.id[3:]][3]) + " studies"
            gene_desc.molecule_expression_cluster_description = expression_sent_enriched + postfix
            joined_sent.append(expression_sent_affected)


def set_do_sentence_and_set_stats(df: DataManager, conf_parser, do_sent_gen_common_props,
                                  do_via_orth_sent_gen_common_props, do_sent_common_props,
                                  do_via_orth_sent_common_props,
                                  gene, gene_desc, joined_sent):
    # This module contains two sub modules for disease and disease via orthology data, each of which is treated as
    # a separate module with different parameters
    do_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                 priority_list=conf_parser.get_do_annotations_priority())
    do_sentence_generator = SentenceGenerator(
        annotations=do_annotations, ontology=df.do_ontology,
        **do_sent_gen_common_props)

    raw_disease_sent = do_sentence_generator.get_sentences(aspect='D', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=False,
                                                           **do_sent_common_props)
    disease_sent = "; ".join([sentence.text for sentence in raw_disease_sent])
    # Disease via orthology module
    do_via_orth_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                          priority_list=conf_parser.
                                                          get_do_via_orth_annotations_priority())
    do_via_orth_sentence_generator = SentenceGenerator(
        annotations=do_via_orth_annotations, ontology=df.do_ontology,
        **do_via_orth_sent_gen_common_props)

    raw_disease_via_orth_sent = do_via_orth_sentence_generator.get_sentences(aspect='D',
                                                                             merge_groups_with_same_prefix=True,
                                                                             keep_only_best_group=False,
                                                                             **do_via_orth_sent_common_props)
    disease_via_orth_sent = "; ".join([sentence.text for sentence in raw_disease_via_orth_sent])
    dis_sent_arr = []
    if len(disease_sent) > 0:
        dis_sent_arr.append(disease_sent)
    if len(disease_via_orth_sent) > 0:
        dis_sent_arr.append(disease_via_orth_sent)
    complete_disease_sent = "; ".join(dis_sent_arr)
    if complete_disease_sent and len(complete_disease_sent) > 0:
        gene_desc.do_description = complete_disease_sent[0].upper() + complete_disease_sent[1:]
        joined_sent.append(complete_disease_sent)
    else:
        gene_desc.do_description = None
    gene_desc.stats.total_number_do_annotations = len(do_annotations) + len(do_via_orth_annotations)
    gene_desc.stats.set_initial_do_ids = [term_id for terms in do_sentence_generator.terms_groups.values() for
                                          tvalues in terms.values() for term_id in tvalues]
    gene_desc.stats.set_initial_do_ids.extend([term_id for terms in
                                               do_via_orth_sentence_generator.terms_groups.values() for
                                               tvalues in terms.values() for term_id in tvalues])
    gene_desc.stats.set_final_do_ids = [term_id for sentence in raw_disease_sent for term_id in
                                        sentence.terms_ids]
    gene_desc.stats.set_final_do_ids.extend([term_id for sentence in raw_disease_sent for term_id in
                                             sentence.terms_ids])
    experimental_disease_sent = "; ".join([sentence.text for sentence in raw_disease_sent if
                                           sentence.evidence_group == "EXPERIMENTAL"])
    if experimental_disease_sent:
        gene_desc.do_experimental_description = experimental_disease_sent
    biomarker_disease_sent = "; ".join([sentence.text for sentence in raw_disease_sent if
                                        sentence.evidence_group == "BIOMARKER"])
    if biomarker_disease_sent:
        gene_desc.do_biomarker_description = biomarker_disease_sent
    orthology_disease_sent = "; ".join([sentence.text for sentence in raw_disease_via_orth_sent if
                                        sentence.evidence_group == "ORTHOLOGY_BASED"])
    if orthology_disease_sent:
        gene_desc.do_orthology_description = orthology_disease_sent
    if "(multiple)" in complete_disease_sent:
        gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms = \
            complete_disease_sent.count("(multiple)")


def set_information_poor_sentence(orthologs_sp_fullname, selected_orthologs, ensembl_hgnc_ids_map, conf_parser,
                                  human_df_agr, go_sent_gen_common_props, go_sent_common_props, gene_desc, df, gene,
                                  joined_sent):
    human_func_sent = None
    if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
        # human_orthologs = df.get_all_orthologs_for_gene(gene_id=gene.id, organism="Homo sapiens")
        best_orth = get_best_human_ortholog_for_info_poor(selected_orthologs, ensembl_hgnc_ids_map,
                                                          conf_parser.get_go_annotations_priority(), human_df_agr,
                                                          go_sent_gen_common_props)
        if best_orth:
            best_orth = "RGD:" + best_orth
            human_go_annotations = human_df_agr.get_annotations_for_gene(
                gene_id=best_orth, annot_type=DataType.GO, priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                                                          "HTP", "HDA", "HMP", "HGI", "HEP"))
            human_go_sent_generator = SentenceGenerator(annotations=human_go_annotations,
                                                        ontology=human_df_agr.go_ontology,
                                                        **go_sent_gen_common_props)
            raw_human_func_sent = human_go_sent_generator.get_sentences(aspect='F',
                                                                        merge_groups_with_same_prefix=True,
                                                                        keep_only_best_group=True,
                                                                        **go_sent_common_props)
            human_func_sent = " and ".join([sentence.text for sentence in raw_human_func_sent])
            if human_func_sent:
                gene_desc.human_gene_function_description = "human " + \
                                                            human_df_agr.go_associations.subject_label_map[
                                                                best_orth] + " " + human_func_sent
                joined_sent.append("human " + human_df_agr.go_associations.subject_label_map[best_orth] + " " +
                                   human_func_sent)
    if not human_func_sent:
        protein_domains = df.protein_domains[gene.id[3:]]
        if protein_domains:
            dom_word = "domain"
            if len(protein_domains) > 1:
                dom_word = "domains"
            joined_sent.append("is predicted to encode a protein with the following " + dom_word + ": " +
                               concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for
                                                                    ptdom in protein_domains]))
            gene_desc.protein_domain_description = \
                "is predicted to encode a protein with the following " + dom_word + ": " + \
                concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for ptdom in
                                                     protein_domains])


def set_sister_species_sentence(df, sister_sp_fullname, sister_df, go_sent_gen_common_props, go_sent_common_props,
                                species, organism, gene_desc, joined_sent, gene):
    best_ortholog = df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
        ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                             "HEP"])[0][0]
    sister_sentences_generator = SentenceGenerator(sister_df.get_annotations_for_gene(
        annot_type=DataType.GO, gene_id="WB:" + best_ortholog[0],
        priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")),
        ontology=df.go_ontology, **go_sent_gen_common_props)
    sister_proc_sent = " and ".join([sentence.text for sentence in sister_sentences_generator.get_sentences(
        aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)])
    if sister_proc_sent:
        gene_desc.sister_species_description = "in " + species[species[organism]["main_sister_species"]]["name"] + \
                                               ", " + best_ortholog[1] + " " + sister_proc_sent
        joined_sent.append("in " + species[species[organism]["main_sister_species"]]["name"] + ", " +
                           best_ortholog[1] + " " + sister_proc_sent)


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.yml", help="configuration file. Default ./config.yaml")
    parser.add_argument("-C", "--use-cache", dest="use_cache", action="store_true", default=False,
                        help="Use cached source files from cache_location specified in config file. Download them from "
                             "raw_file_source (configured in config file) if not yet cached")
    parser.add_argument("-l", "--log-file", metavar="log_file", dest="log_file", type=str, default=None,
                        help="path to the log file to generate. Default ./genedescriptions.log")
    parser.add_argument("-L", "--log-level", dest="log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], help="set the logging level")
    parser.add_argument("-o", "--output-formats", metavar="output_formats", dest="output_formats", type=str, nargs="+",
                        default=["ace", "txt", "json", "tsv"], help="file formats to generate. Accepted values "
                                                                    "are: ace, txt, json, tsv")

    args = parser.parse_args()
    conf_parser = GenedescConfigParser(args.config_file)
    logging.basicConfig(filename=args.log_file, level=args.log_level, format='%(asctime)s - %(name)s - %(levelname)s:'
                                                                             '%(message)s')

    logger = logging.getLogger("WB Gene Description Pipeline")

    go_sent_gen_common_props = {"evidence_groups_priority_list": conf_parser.get_go_evidence_groups_priority_list(),
                                "prepostfix_sentences_map": conf_parser.get_go_prepostfix_sentences_map(),
                                "prepostfix_special_cases_sent_map":
                                    conf_parser.get_go_prepostfix_special_cases_sent_map(),
                                "evidence_codes_groups_map": conf_parser.get_go_evidence_codes_groups_map()}
    go_sent_common_props = {"remove_parent_terms": conf_parser.get_go_remove_parents_if_children_are_present(),
                            "remove_child_terms": conf_parser.get_go_remove_children_if_parent_is_present(),
                            "merge_num_terms_threshold": conf_parser.get_go_trim_min_num_terms(),
                            "merge_max_num_terms": conf_parser.get_go_max_num_terms(),
                            "merge_min_distance_from_root": conf_parser.get_go_trim_min_distance_from_root(),
                            "truncate_others_generic_word": conf_parser.get_go_truncate_others_aggregation_word(),
                            "truncate_others_aspect_words": conf_parser.get_go_truncate_others_terms(),
                            "add_multiple_if_covers_more_children": False,
                            "blacklisted_ancestors": conf_parser.get_go_terms_exclusion_list()}
    do_sent_gen_common_props = {"evidence_groups_priority_list": conf_parser.get_do_evidence_groups_priority_list(),
                                "prepostfix_sentences_map": conf_parser.get_do_prepostfix_sentences_map(),
                                "prepostfix_special_cases_sent_map": None,
                                "evidence_codes_groups_map": conf_parser.get_do_evidence_codes_groups_map()}
    do_via_orth_sent_gen_common_props = {
        "evidence_groups_priority_list": conf_parser.get_do_via_orth_evidence_groups_priority_list(),
        "prepostfix_sentences_map": conf_parser.get_do_via_orth_prepostfix_sentences_map(),
        "prepostfix_special_cases_sent_map": None,
        "evidence_codes_groups_map": conf_parser.get_do_via_orth_evidence_codes_groups_map()}
    do_sent_common_props = {"remove_parent_terms": conf_parser.get_do_remove_parents_if_children_are_present(),
                            "remove_child_terms": conf_parser.get_do_remove_children_if_parent_is_present(),
                            "merge_num_terms_threshold": conf_parser.get_do_trim_min_num_terms(),
                            "merge_max_num_terms": conf_parser.get_do_max_num_terms(),
                            "merge_min_distance_from_root": conf_parser.get_do_trim_min_distance_from_root(),
                            "truncate_others_generic_word": conf_parser.get_do_truncate_others_aggregation_word(),
                            "truncate_others_aspect_words": conf_parser.get_do_truncate_others_terms(),
                            "add_multiple_if_covers_more_children": True,
                            "blacklisted_ancestors": conf_parser.get_do_terms_exclusion_list()}
    do_via_orth_sent_common_props = {"remove_parent_terms": conf_parser.get_do_remove_parents_if_children_are_present(),
                                     "remove_child_terms": conf_parser.get_do_remove_children_if_parent_is_present(),
                                     "merge_num_terms_threshold": conf_parser.get_do_trim_min_num_terms(),
                                     "merge_max_num_terms": conf_parser.get_do_max_num_terms(),
                                     "merge_min_distance_from_root": conf_parser.get_do_trim_min_distance_from_root(),
                                     "truncate_others_generic_word": conf_parser.get_do_truncate_others_aggregation_word(),
                                     "truncate_others_aspect_words": conf_parser.get_do_truncate_others_terms(),
                                     "add_multiple_if_covers_more_children": True,
                                     "blacklisted_ancestors": conf_parser.get_do_via_orth_terms_exclusion_list()}
    expr_sent_gen_common_props = {
        "evidence_groups_priority_list": conf_parser.get_expression_evidence_groups_priority_list(),
        "prepostfix_sentences_map": conf_parser.get_expression_prepostfix_sentences_map(),
        "prepostfix_special_cases_sent_map": None,
        "evidence_codes_groups_map": conf_parser.get_expression_evidence_codes_groups_map()}
    expr_sent_common_props = {
        "remove_parent_terms": conf_parser.get_expression_remove_parents_if_children_are_present(),
        "remove_child_terms": conf_parser.get_expression_remove_children_if_parent_is_present(),
        "merge_num_terms_threshold": conf_parser.get_expression_trim_min_num_terms(),
        "merge_max_num_terms": conf_parser.get_expression_max_num_terms(),
        "merge_min_distance_from_root": conf_parser.get_expression_trim_min_distance_from_root(),
        "truncate_others_generic_word": conf_parser.get_expression_truncate_others_aggregation_word(),
        "truncate_others_aspect_words": conf_parser.get_expression_truncate_others_terms(),
        "add_multiple_if_covers_more_children": False, "rename_cell": True,
        "blacklisted_ancestors": conf_parser.get_expression_terms_exclusion_list()}

    organisms_list = conf_parser.get_wb_organisms_to_process()
    human_genes_props = DataManager.get_human_gene_props()
    ensembl_hgnc_ids_map = DataManager.get_ensembl_hgnc_ids_map()
    for organism in organisms_list:
        logger.info("Processing organism " + organism)
        species = conf_parser.get_wb_species()
        df, sister_df, df_agr, orthologs_sp_fullname, sister_sp_fullname = load_data(
            species=species, organism=organism, conf_parser=conf_parser)
        desc_writer = DescriptionsWriter()
        desc_writer.overall_properties.species = organism
        desc_writer.overall_properties.release_version = conf_parser.get_release("wb_data_fetcher")[0:-1] + str(
            int(conf_parser.get_release("wb_data_fetcher")[-1]) + 1)
        desc_writer.overall_properties.date = datetime.date.today().strftime("%B %d, %Y")
        for gene in df.get_gene_data():
            logger.debug("Generating description for gene " + gene.name)

            gene_desc = GeneDesc(gene_id=gene.id, gene_name=gene.name,
                                 publications=", ".join([annot["publication"] for annot in df.get_annotations_for_gene(
                                     gene.id, annot_type=DataType.GO,
                                     priority_list=conf_parser.get_go_evidence_groups_priority_list())]),
                                 refs=", ".join([annot["refs"] for annot in df.get_annotations_for_gene(
                                     gene.id, annot_type=DataType.GO,
                                     priority_list=conf_parser.get_go_evidence_groups_priority_list())]))
            joined_sent = []

            selected_orthologs = add_orthology_sentence(df=df, gene=gene, orthologs_sp_fullname=orthologs_sp_fullname,
                                                        human_genes_props=human_genes_props, gene_desc=gene_desc,
                                                        conf_parser=conf_parser, joined_sent=joined_sent)
            go_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                         priority_list=conf_parser.get_go_annotations_priority())
            go_sent_gen_common_props_exp = go_sent_gen_common_props.copy()
            go_sent_gen_common_props_exp["evidence_codes_groups_map"] = {
                evcode: group for evcode, group in go_sent_gen_common_props["evidence_codes_groups_map"].items() if
                "EXPERIMENTAL" in go_sent_gen_common_props["evidence_codes_groups_map"][evcode]}
            go_sent_generator_exp = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                                      **go_sent_gen_common_props_exp)
            go_sent_generator = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                                  **go_sent_gen_common_props)
            set_initial_go_stats(go_annotations=go_annotations, go_sent_generator=go_sent_generator,
                                 go_sent_generator_exp=go_sent_generator_exp, gene_desc=gene_desc,
                                 conf_parser=conf_parser)
            set_go_function_sentence_and_set_stats(go_sent_generator=go_sent_generator,
                                                   go_sent_generator_exp=go_sent_generator_exp,
                                                   go_sent_common_props=go_sent_common_props, gene_desc=gene_desc,
                                                   joined_sent=joined_sent)
            set_go_process_sentence_and_set_stats(go_sent_generator=go_sent_generator,
                                                  go_sent_common_props=go_sent_common_props, gene_desc=gene_desc,
                                                  joined_sent=joined_sent)
            set_go_component_sentence_and_set_stats(go_sent_generator=go_sent_generator,
                                                    go_sent_generator_exp=go_sent_generator_exp,
                                                    go_sent_common_props=go_sent_common_props, gene_desc=gene_desc,
                                                    joined_sent=joined_sent)
            expr_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.EXPR,
                                                           priority_list=conf_parser.get_expression_annotations_priority())
            expr_sentence_generator = SentenceGenerator(annotations=expr_annotations, ontology=df.expression_ontology,
                                                        **expr_sent_gen_common_props)
            gene_desc.stats.set_initial_expression_ids = [elem for key, sets in
                                                          expr_sentence_generator.terms_groups[
                                                              ('A', 'Verified')].items() for
                                                          elem in sets if ('A', key, 'Verified') in
                                                          conf_parser.get_expression_prepostfix_sentences_map()]
            set_expression_sentence(expr_sentence_generator=expr_sentence_generator,
                                    expr_sent_common_props=expr_sent_common_props, gene_desc=gene_desc, gene=gene,
                                    joined_sent=joined_sent, df=df)

            set_do_sentence_and_set_stats(df=df, conf_parser=conf_parser, do_sent_common_props=do_sent_common_props,
                                          do_via_orth_sent_common_props=do_via_orth_sent_common_props, gene=gene,
                                          gene_desc=gene_desc, joined_sent=joined_sent,
                                          do_sent_gen_common_props=do_sent_gen_common_props,
                                          do_via_orth_sent_gen_common_props=do_via_orth_sent_gen_common_props)
            if not gene_desc.go_description:
                set_information_poor_sentence(orthologs_sp_fullname=orthologs_sp_fullname,
                                              selected_orthologs=selected_orthologs,
                                              ensembl_hgnc_ids_map=ensembl_hgnc_ids_map,
                                              conf_parser=conf_parser, human_df_agr=df_agr,
                                              go_sent_gen_common_props=go_sent_gen_common_props,
                                              go_sent_common_props=go_sent_common_props, gene_desc=gene_desc, df=df,
                                              gene=gene,
                                              joined_sent=joined_sent)
            if conf_parser.get_data_fetcher() == "wb_data_fetcher" and "main_sister_species" in species[organism] and \
                    species[organism]["main_sister_species"] and df.get_best_orthologs_for_gene(
                gene.id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
                ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                                     "HEP"])[0]:
                set_sister_species_sentence(df=df, sister_sp_fullname=sister_sp_fullname, sister_df=sister_df,
                                            go_sent_gen_common_props=go_sent_gen_common_props,
                                            go_sent_common_props=go_sent_common_props,
                                            species=species, organism=organism, gene_desc=gene_desc,
                                            joined_sent=joined_sent,
                                            gene=gene)

            if len(joined_sent) > 0:
                desc = "; ".join(joined_sent) + "."
                if len(desc) > 0:
                    gene_desc.description = gene.name + " " + desc
            else:
                gene_desc.description = None
            desc_writer.add_gene_desc(gene_desc)
        if "json" in args.output_formats:
            desc_writer.write_json(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                                organism + ".json"), pretty=True, include_single_gene_stats=True)
        if "txt" in args.output_formats:
            desc_writer.write_plain_text(os.path.join(
                conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()), organism + ".txt"))
        if "tsv" in args.output_formats:
            desc_writer.write_tsv(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                               organism + ".tsv"))
        if "ace" in args.output_formats:
            curators = ["WBPerson324", "WBPerson37462"]
            release_version = conf_parser.get_release("wb_data_fetcher")
            desc_writer.write_ace(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                               organism + ".ace"), curators, release_version)


if __name__ == '__main__':
    main()
