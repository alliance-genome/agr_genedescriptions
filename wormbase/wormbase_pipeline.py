#!/usr/bin/env python3

import argparse
import datetime
import logging
import os

from genedescriptions.commons import DataType, Module
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import WBDataManager, DataManager, ExpressionClusterType, ExpressionClusterFeature
from genedescriptions.descriptions_and_stats import GeneDescription
from genedescriptions.descriptions_generator import SentenceGenerator
from genedescriptions.descriptions_writer import DescriptionsWriter
from genedescriptions.sentence_generation_functions import generate_ortholog_sentence_wormbase_human, \
    generate_ortholog_sentence_wormbase_non_c_elegans, concatenate_words_with_oxford_comma, \
    get_best_human_ortholog_for_info_poor


def load_data(species, organism, conf_parser: GenedescConfigParser):
    logger = logging.getLogger("WB Gene Description Pipeline - Data loader")
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
    ec_anatomy_prefix = species[organism]["ec_anatomy_prefix"] if "ec_anatomy_prefix" in species[organism] else None
    ec_molreg_prefix = species[organism]["ec_molreg_prefix"] if "ec_molreg_prefix" in species[organism] else None
    df = WBDataManager(raw_files_source=conf_parser.get_wb_raw_file_sources(),
                       release_version=conf_parser.get_wb_release(),
                       species=organism, project_id=species[organism]["project_id"],
                       cache_location=conf_parser.get_cache_dir(), do_relations=None,
                       go_relations=["subClassOf", "BFO:0000050"], sister_sp_fullname=sister_sp_fullname,
                       expression_cluster_anatomy_prefix=ec_anatomy_prefix,
                       expression_cluster_molreg_prefix=ec_molreg_prefix)
    if organism == "c_elegans":
        df_agr = DataManager(go_relations=["subClassOf", "BFO:0000050"], do_relations=None)
        df_agr.load_ontology_from_file(ontology_type=DataType.GO,
                                       ontology_url=conf_parser.get_wb_human_orthologs_go_ontology(),
                                       ontology_cache_path=os.path.join(
                                           conf_parser.get_cache_dir(),
                                           "wormbase_agr_human", "go_ontology.obo"),
                                       terms_replacement_regex=conf_parser.get_module_simple_property(
                                           module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS))
        df_agr.load_associations_from_file(associations_type=DataType.GO,
                                           associations_url=conf_parser.get_wb_human_orthologs_go_associations(),
                                           associations_cache_path=os.path.join(
                                               conf_parser.get_cache_dir(), "wormbase_agr_human", "go_assoc.daf.gz"),
                                           exclusion_list=conf_parser.get_module_simple_property(
                                               module=Module.GO, prop=ConfigModuleProperty.EXCLUDE_TERMS))
    if "main_sister_species" in species[organism] and species[organism]["main_sister_species"]:
        sister_df = WBDataManager(raw_files_source=conf_parser.get_wb_raw_file_sources(),
                                  release_version=conf_parser.get_wb_release(),
                                  species=species[organism]["main_sister_species"],
                                  project_id=species[species[organism]["main_sister_species"]]["project_id"],
                                  cache_location=conf_parser.get_cache_dir(), do_relations=None,
                                  go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading all data for sister species")
        sister_df.load_all_data_from_file(
            go_terms_replacement_regex=conf_parser.get_module_simple_property(
                module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS),
            go_terms_exclusion_list=conf_parser.get_module_simple_property(
                module=Module.GO, prop=ConfigModuleProperty.EXCLUDE_TERMS),
            do_terms_replacement_regex=conf_parser.get_module_simple_property(
                module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.RENAME_TERMS),
            do_terms_exclusion_list=conf_parser.get_module_simple_property(
                module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.EXCLUDE_TERMS))
    logger.info("Loading all data for main species")
    df.load_all_data_from_file(
        go_terms_replacement_regex=conf_parser.get_module_simple_property(
            module=Module.GO, prop=ConfigModuleProperty.RENAME_TERMS),
        go_terms_exclusion_list=conf_parser.get_module_simple_property(
            module=Module.GO, prop=ConfigModuleProperty.EXCLUDE_TERMS),
        do_terms_replacement_regex=conf_parser.get_module_simple_property(
            module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.RENAME_TERMS),
        do_terms_exclusion_list=conf_parser.get_module_simple_property(
            module=Module.DO_EXP_AND_BIO, prop=ConfigModuleProperty.EXCLUDE_TERMS))
    if organism == "c_elegans":
        df.load_ontology_from_file(ontology_type=DataType.EXPR, ontology_url=df.expression_ontology_url,
                                   ontology_cache_path=df.expression_ontology_cache_path,
                                   terms_replacement_regex=conf_parser.get_module_simple_property(
                                       module=Module.EXPRESSION, prop=ConfigModuleProperty.RENAME_TERMS))
        df.load_associations_from_file(associations_type=DataType.EXPR,
                                       associations_url=df.expression_associations_url,
                                       associations_cache_path=df.expression_associations_cache_path,
                                       exclusion_list=conf_parser.get_module_simple_property(
                                           module=Module.EXPRESSION, prop=ConfigModuleProperty.EXCLUDE_TERMS))
    df.load_expression_cluster_data()
    return df, sister_df, df_agr, orthologs_sp_fullname, sister_sp_fullname


def set_orthology_sentence(df, gene, orthologs_sp_fullname, human_genes_props, gene_desc, tpc_token):
    best_orthologs, selected_orth_name = df.get_best_orthologs_for_gene(gene.id,
                                                                        orth_species_full_name=orthologs_sp_fullname)
    selected_orthologs = []
    if best_orthologs:
        gene_desc.stats.set_best_orthologs = [orth[0] for orth in best_orthologs]
        if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
            sel_orthologs, orth_sent = generate_ortholog_sentence_wormbase_human(best_orthologs, human_genes_props)
            selected_orthologs = [orth for orth in best_orthologs if orth[1] in sel_orthologs]
        else:
            orth_sent = generate_ortholog_sentence_wormbase_non_c_elegans(best_orthologs, selected_orth_name, tpc_token)
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.ORTHOLOGY, description=orth_sent)
    return selected_orthologs


def set_go_sentences(df, go_sent_gen_common_props, go_sent_gen_common_props_exp, go_sent_common_props, gene,
                     conf_parser: GenedescConfigParser, gene_desc):
    go_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                 priority_list=conf_parser.get_annotations_priority(module=Module.GO))

    go_sent_generator_exp = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                              **go_sent_gen_common_props_exp)
    go_sent_generator = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                          **go_sent_gen_common_props)
    contributes_to_module_sentences = go_sent_generator.get_module_sentences(
        aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True, keep_only_best_group=True,
        **go_sent_common_props)
    if contributes_to_module_sentences.contains_sentences():
        func_module_sentences = go_sent_generator_exp.get_module_sentences(
            aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=func_module_sentences, module=Module.GO_FUNCTION, annotations=go_annotations)
    else:
        func_module_sentences = go_sent_generator.get_module_sentences(
            aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=func_module_sentences, module=Module.GO_FUNCTION, annotations=go_annotations)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=contributes_to_module_sentences, module=Module.GO_FUNCTION, annotations=go_annotations)
    proc_module_sentences = go_sent_generator.get_module_sentences(aspect='P', merge_groups_with_same_prefix=True,
                                                                   keep_only_best_group=True, **go_sent_common_props)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences, module=Module.GO_PROCESS, annotations=go_annotations)
    colocalizes_with_module_sentences = go_sent_generator.get_module_sentences(
        aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True, **go_sent_common_props)
    if colocalizes_with_module_sentences.contains_sentences():
        comp_module_sentences = go_sent_generator_exp.get_module_sentences(
            aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT, annotations=go_annotations)
    else:
        comp_module_sentences = go_sent_generator.get_module_sentences(
            aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT, annotations=go_annotations)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=colocalizes_with_module_sentences,
                                                               module=Module.GO_COMPONENT, annotations=go_annotations)
    gene_desc.set_initial_stats(module=Module.GO_FUNCTION, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_PROCESS, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_COMPONENT, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)


def set_expression_sentence(expr_sent_gen_common_props, expr_sent_common_props, gene_desc, gene, df,
                            conf_parser: GenedescConfigParser):
    expr_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.EXPR,
                                                   priority_list=conf_parser.get_annotations_priority(
                                                       module=Module.EXPRESSION))
    expr_sentence_generator = SentenceGenerator(annotations=expr_annotations, ontology=df.expression_ontology,
                                                **expr_sent_gen_common_props)
    expression_module_sentences = expr_sentence_generator.get_module_sentences(
        aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True, keep_only_best_group=False,
        **expr_sent_common_props)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=expression_module_sentences,
                                                               module=Module.EXPRESSION)
    gene_desc.set_initial_stats(module=Module.EXPRESSION, sentence_generator=expr_sentence_generator)
    # Information poor genes
    if not gene_desc.description:
        ec_gene_id = gene.id[3:]
        ec_anatomy_studies = df.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                               expression_cluster_type=ExpressionClusterType.ANATOMY,
                                                               feature=ExpressionClusterFeature.STUDIES)
        ec_anatomy_terms = df.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                             feature=ExpressionClusterFeature.TERMS,
                                                             expression_cluster_type=ExpressionClusterType.ANATOMY)
        if df.expression_ontology is not None:
            expression_enriched_module_sentences = expr_sentence_generator.get_module_sentences(
                aspect='A', qualifier="Enriched", merge_groups_with_same_prefix=True, keep_only_best_group=False,
                **expr_sent_common_props)
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.EXPRESSION_CLUSTER_GENE,
                description=expression_enriched_module_sentences.get_description(),
                additional_postfix_terms_list=ec_anatomy_studies, additional_postfix_final_word="studies",
                use_single_form=True)
        elif ec_anatomy_terms:
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.EXPRESSION_CLUSTER_ANATOMY,
                description="is enriched in " + concatenate_words_with_oxford_comma(ec_anatomy_terms) + " based on ",
                additional_postfix_terms_list=ec_anatomy_studies,
                additional_postfix_final_word="studies", use_single_form=True)
        ec_molreg_terms = df.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                            expression_cluster_type=ExpressionClusterType.MOLREG,
                                                            feature=ExpressionClusterFeature.TERMS)
        ec_molereg_studies = df.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                               feature=ExpressionClusterFeature.STUDIES,
                                                               expression_cluster_type=ExpressionClusterType.MOLREG)

        if df.expression_ontology is None and ec_molreg_terms:
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.EXPRESSION_CLUSTER_MOLECULE,
                description="is affected by " + concatenate_words_with_oxford_comma(ec_molreg_terms) + " based on ",
                additional_postfix_terms_list=ec_molereg_studies,
                additional_postfix_final_word="studies", use_single_form=True)


def set_do_sentence(df: DataManager, conf_parser: GenedescConfigParser, do_sent_gen_common_props,
                    do_via_orth_sent_gen_common_props, do_sent_common_props, do_via_orth_sent_common_props, gene,
                    gene_desc):
    do_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                 priority_list=conf_parser.get_annotations_priority(
                                                     module=Module.DO_EXP_AND_BIO))
    # Experimental
    do_sent_gen_common_props_exp = do_sent_gen_common_props.copy()
    do_sent_gen_common_props_exp["evidence_codes_groups_map"] = {
        evcode: group for evcode, group in do_sent_gen_common_props["evidence_codes_groups_map"].items() if
        "EXPERIMENTAL" in do_sent_gen_common_props["evidence_codes_groups_map"][evcode]}
    do_sentence_exp_generator = SentenceGenerator(annotations=do_annotations, ontology=df.do_ontology,
                                                  **do_sent_gen_common_props_exp)
    disease_exp_module_sentences = do_sentence_exp_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, **do_sent_common_props)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_EXPERIMENTAL,
                                                               module_sentences=disease_exp_module_sentences,
                                                               annotations=do_annotations)
    # Biomarker
    do_sent_gen_common_props_bio = do_sent_gen_common_props.copy()
    do_sent_gen_common_props_bio["evidence_codes_groups_map"] = {
        evcode: group for evcode, group in do_sent_gen_common_props["evidence_codes_groups_map"].items() if
        "BIOMARKER" in do_sent_gen_common_props["evidence_codes_groups_map"][evcode]}
    do_sentence_bio_generator = SentenceGenerator(annotations=do_annotations, ontology=df.do_ontology,
                                                  **do_sent_gen_common_props_bio)
    disease_bio_module_sentences = do_sentence_bio_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, **do_sent_common_props)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_BIOMARKER,
                                                               module_sentences=disease_bio_module_sentences,
                                                               annotations=do_annotations)
    # Disease via orthology module
    do_via_orth_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                          priority_list=conf_parser.
                                                          get_annotations_priority(module=Module.DO_ORTHOLOGY))
    do_via_orth_sentence_generator = SentenceGenerator(annotations=do_via_orth_annotations, ontology=df.do_ontology,
                                                       **do_via_orth_sent_gen_common_props)
    disease_via_orth_module_sentences = do_via_orth_sentence_generator.get_module_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, **do_via_orth_sent_common_props)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_ORTHOLOGY,
                                                               module_sentences=disease_via_orth_module_sentences,
                                                               annotations=do_via_orth_annotations)
    gene_desc.set_initial_stats(module=Module.DO_EXPERIMENTAL, sentence_generator=do_sentence_exp_generator)
    gene_desc.set_initial_stats(module=Module.DO_BIOMARKER, sentence_generator=do_sentence_bio_generator)
    gene_desc.set_initial_stats(module=Module.DO_ORTHOLOGY, sentence_generator=do_via_orth_sentence_generator)


def set_information_poor_sentence(orthologs_sp_fullname, selected_orthologs, ensembl_hgnc_ids_map,
                                  conf_parser: GenedescConfigParser, human_df_agr, go_sent_gen_common_props,
                                  go_sent_common_props, gene_desc, df, gene):
    human_func_sent = None
    if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
        best_orth = get_best_human_ortholog_for_info_poor(selected_orthologs, ensembl_hgnc_ids_map,
                                                          conf_parser.get_annotations_priority(module=Module.GO),
                                                          human_df_agr, go_sent_gen_common_props)
        if best_orth:
            best_orth = "RGD:" + best_orth
            human_go_annotations = human_df_agr.get_annotations_for_gene(
                gene_id=best_orth, annot_type=DataType.GO, priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                                                          "HTP", "HDA", "HMP", "HGI", "HEP"))
            human_go_sent_generator = SentenceGenerator(annotations=human_go_annotations,
                                                        ontology=human_df_agr.go_ontology,
                                                        **go_sent_gen_common_props)
            human_func_module_sentences = human_go_sent_generator.get_module_sentences(
                aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
            human_func_sent = human_func_module_sentences.get_description()
            if human_func_sent:
                gene_desc.set_or_extend_module_description_and_final_stats(
                    module=Module.INFO_POOR_HUMAN_FUNCTION, description="human " +
                                                                        human_df_agr.go_associations.subject_label_map[
                                                                            best_orth] + " " + human_func_sent)
    if not human_func_sent:
        protein_domains = df.protein_domains[gene.id[3:]]
        if protein_domains:
            dom_word = "domain"
            if len(protein_domains) > 1:
                dom_word = "domains"
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.INFO_POOR_HUMAN_FUNCTION,
                description="is predicted to encode a protein with the following " + dom_word + ": " +
                            concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for
                                                                 ptdom in protein_domains]))


def set_sister_species_sentence(df, sister_sp_fullname, sister_df, go_sent_gen_common_props, go_sent_common_props,
                                species, organism, gene_desc, gene):
    best_ortholog = df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
        ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                             "HEP"])[0][0]
    sister_sentences_generator = SentenceGenerator(sister_df.get_annotations_for_gene(
        annot_type=DataType.GO, gene_id="WB:" + best_ortholog[0],
        priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")),
        ontology=df.go_ontology, **go_sent_gen_common_props)
    sister_sp_module_sentences = sister_sentences_generator.get_module_sentences(
        aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)\
        .get_description()
    gene_desc.set_or_extend_module_description_and_final_stats(
        module=Module.SISTER_SP, description="in " + species[species[organism]["main_sister_species"]]["name"] +
                                             ", " + best_ortholog[1] + " " + sister_sp_module_sentences)


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
    parser.add_argument("-t", "--textpressoapi-token", metavar="textpresso_token", dest="textpresso_token", type=str,
                        help="Texpresso api token")
    parser.add_argument("-o", "--output-formats", metavar="output_formats", dest="output_formats", type=str, nargs="+",
                        default=["ace", "txt", "json", "tsv"], help="file formats to generate. Accepted values "
                                                                    "are: ace, txt, json, tsv")

    args = parser.parse_args()
    conf_parser = GenedescConfigParser(args.config_file)
    logging.basicConfig(filename=args.log_file, level=args.log_level, format='%(asctime)s - %(name)s - %(levelname)s:'
                                                                             '%(message)s')

    logger = logging.getLogger("WB Gene Description Pipeline")
    organisms_list = conf_parser.get_wb_organisms_to_process()
    human_genes_props = DataManager.get_human_gene_props()
    ensembl_hgnc_ids_map = DataManager.get_ensembl_hgnc_ids_map()
    go_sent_gen_common_props = conf_parser.get_sentence_generator_common_properties(Module.GO)
    go_sent_gen_common_props_exp = go_sent_gen_common_props.copy()
    go_sent_common_props = conf_parser.get_sentence_common_properties(module=Module.GO)
    expr_sent_gen_common_props = conf_parser.get_sentence_generator_common_properties(module=Module.EXPRESSION)
    expr_sent_common_props = conf_parser.get_sentence_common_properties(module=Module.EXPRESSION)
    do_sent_gen_common_props = conf_parser.get_sentence_generator_common_properties(module=Module.DO_EXP_AND_BIO)
    do_sent_common_props = conf_parser.get_sentence_common_properties(module=Module.DO_EXP_AND_BIO)
    do_via_orth_sent_gen_common_props = conf_parser.get_sentence_generator_common_properties(module=Module.DO_ORTHOLOGY)
    do_via_orth_sent_common_props = conf_parser.get_sentence_common_properties(module=Module.DO_ORTHOLOGY)
    go_sent_gen_common_props_exp["evidence_codes_groups_map"] = {
        evcode: group for evcode, group in go_sent_gen_common_props["evidence_codes_groups_map"].items() if
        "EXPERIMENTAL" in go_sent_gen_common_props["evidence_codes_groups_map"][evcode]}
    for organism in organisms_list:
        logger.info("Processing organism " + organism)
        species = conf_parser.get_wb_organisms_info()
        df, sister_df, df_agr, orthologs_sp_fullname, sister_sp_fullname = load_data(
            species=species, organism=organism, conf_parser=conf_parser)
        desc_writer = DescriptionsWriter()
        desc_writer.overall_properties.species = organism
        desc_writer.overall_properties.release_version = conf_parser.get_wb_release()[0:-1] + str(
            int(conf_parser.get_wb_release()[-1]) + 1)
        desc_writer.overall_properties.date = datetime.date.today().strftime("%B %d, %Y")
        for gene in df.get_gene_data():
            logger.debug("Generating description for gene " + gene.name)

            gene_desc = GeneDescription(gene_id=gene.id, gene_name=gene.name, add_gene_name=True)
            selected_orthologs = set_orthology_sentence(df=df, gene=gene, orthologs_sp_fullname=orthologs_sp_fullname,
                                                        human_genes_props=human_genes_props, gene_desc=gene_desc,
                                                        tpc_token=args.textpresso_token)
            set_go_sentences(df=df, go_sent_gen_common_props=go_sent_gen_common_props,
                             go_sent_gen_common_props_exp=go_sent_gen_common_props_exp,
                             go_sent_common_props=go_sent_common_props, gene=gene, conf_parser=conf_parser,
                             gene_desc=gene_desc)
            set_expression_sentence(expr_sent_gen_common_props=expr_sent_gen_common_props,
                                    expr_sent_common_props=expr_sent_common_props, gene_desc=gene_desc, gene=gene,
                                    df=df, conf_parser=conf_parser)

            set_do_sentence(df=df, conf_parser=conf_parser, do_sent_common_props=do_sent_common_props,
                            do_via_orth_sent_common_props=do_via_orth_sent_common_props, gene=gene,
                            gene_desc=gene_desc, do_sent_gen_common_props=do_sent_gen_common_props,
                            do_via_orth_sent_gen_common_props=do_via_orth_sent_gen_common_props)
            if not gene_desc.go_description:
                set_information_poor_sentence(orthologs_sp_fullname=orthologs_sp_fullname,
                                              selected_orthologs=selected_orthologs,
                                              ensembl_hgnc_ids_map=ensembl_hgnc_ids_map,
                                              conf_parser=conf_parser, human_df_agr=df_agr,
                                              go_sent_gen_common_props=go_sent_gen_common_props,
                                              go_sent_common_props=go_sent_common_props, gene_desc=gene_desc, df=df,
                                              gene=gene)
            if "main_sister_species" in species[organism] and species[organism]["main_sister_species"] and \
                    df.get_best_orthologs_for_gene(gene.id, orth_species_full_name=[sister_sp_fullname],
                                                   sister_species_data_fetcher=sister_df,
                                                   ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP",
                                                                        "HDA", "HMP", "HGI", "HEP"])[0]:
                set_sister_species_sentence(df=df, sister_sp_fullname=sister_sp_fullname, sister_df=sister_df,
                                            go_sent_gen_common_props=go_sent_gen_common_props,
                                            go_sent_common_props=go_sent_common_props,
                                            species=species, organism=organism, gene_desc=gene_desc, gene=gene)
            desc_writer.add_gene_desc(gene_desc)
        if "json" in args.output_formats:
            desc_writer.write_json(os.path.join(conf_parser.get_out_dir(), organism + ".json"),
                                   pretty=True, include_single_gene_stats=True)
        if "txt" in args.output_formats:
            desc_writer.write_plain_text(os.path.join(conf_parser.get_out_dir(), organism + ".txt"))
        if "tsv" in args.output_formats:
            desc_writer.write_tsv(os.path.join(conf_parser.get_out_dir(), organism + ".tsv"))
        if "ace" in args.output_formats:
            curators = ["WBPerson324", "WBPerson37462"]
            release_version = conf_parser.get_wb_release()
            desc_writer.write_ace(os.path.join(conf_parser.get_out_dir(), organism + ".ace"), curators, release_version)


if __name__ == '__main__':
    main()
