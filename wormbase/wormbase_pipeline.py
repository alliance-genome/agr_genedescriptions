#!/usr/bin/env python3

import argparse
import datetime
import logging
import os
from typing import List

from genedescriptions.commons import DataType, Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, ExpressionClusterType, ExpressionClusterFeature
from genedescriptions.gene_description import GeneDescription
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.descriptions_writer import DescriptionsWriter
from genedescriptions.sentence_generation_functions import generate_ortholog_sentence_wormbase_human, \
    generate_ortholog_sentence_wormbase_non_c_elegans, concatenate_words_with_oxford_comma, \
    get_best_human_ortholog_for_info_poor
from wormbase.wb_data_manager import WBDataManager


def load_data(species, organism, conf_parser: GenedescConfigParser):
    logger = logging.getLogger("WB Gene Description Pipeline - Data loader")
    sister_df = None
    df_agr = None
    sister_sp_fullname = ""
    if "main_sister_species" in species[organism] and "full_name" in \
            species[species[organism]["main_sister_species"]]:
        sister_sp_fullname = species[species[organism]["main_sister_species"]]["full_name"]
    orth_fullnames = ""
    if "ortholog" in species[organism] and all(["full_name" in species[ortholog_sp] for ortholog_sp in
                                                species[organism]["ortholog"]]):
        orth_fullnames = [species[ortholog_sp]["full_name"] for ortholog_sp in species[organism]["ortholog"]]
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
    return df, sister_df, df_agr, orth_fullnames, sister_sp_fullname


def set_orthology_sentence(dm: WBDataManager, orth_fullnames: List[str], gene_desc: GeneDescription, tpc_token: str,
                           human_genes_props):
    best_orthologs, selected_orth_name = dm.get_best_orthologs_for_gene(gene_desc.gene_id,
                                                                        orth_species_full_name=orth_fullnames)
    selected_orthologs = []
    if best_orthologs:
        gene_desc.stats.set_best_orthologs = [orth[0] for orth in best_orthologs]
        if len(orth_fullnames) == 1 and orth_fullnames[0] == "Homo sapiens":
            sel_orthologs, orth_sent = generate_ortholog_sentence_wormbase_human(best_orthologs, human_genes_props)
            selected_orthologs = [orth for orth in best_orthologs if orth[1] in sel_orthologs]
        else:
            orth_sent = generate_ortholog_sentence_wormbase_non_c_elegans(best_orthologs, selected_orth_name, tpc_token)
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.ORTHOLOGY, description=orth_sent)
    return selected_orthologs


def set_go_sentences(dm: WBDataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription, gene: Gene):
    go_sent_generator_exp = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.GO,
                                                      module=Module.GO, data_manager=dm,
                                                      config=conf_parser, limit_to_group="EXPERIMENTAL")
    go_sent_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.GO,
                                                  module=Module.GO, data_manager=dm, config=conf_parser)
    contributes_to_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if contributes_to_module_sentences.contains_sentences():
        func_module_sentences = go_sent_generator_exp.get_module_sentences(
            config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=func_module_sentences,
                                                                   module=Module.GO_FUNCTION)
    else:
        func_module_sentences = go_sent_generator.get_module_sentences(
            config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=func_module_sentences, module=Module.GO_FUNCTION)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=contributes_to_module_sentences, module=Module.GO_FUNCTION)
    proc_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences, module=Module.GO_PROCESS)
    colocalizes_with_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if colocalizes_with_module_sentences.contains_sentences():
        comp_module_sentences = go_sent_generator_exp.get_module_sentences(
            config=conf_parser, aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
    else:
        comp_module_sentences = go_sent_generator.get_module_sentences(
            config=conf_parser, aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=colocalizes_with_module_sentences,
                                                               module=Module.GO_COMPONENT)
    gene_desc.set_initial_stats(module=Module.GO_FUNCTION, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_PROCESS, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_COMPONENT, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)


def set_expression_sentence(dm: WBDataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription,
                            gene: Gene):
    expr_sentence_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.GO,
                                                        module=Module.GO, data_manager=dm,
                                                        config=conf_parser, limit_to_group="EXPERIMENTAL")
    expression_module_sentences = expr_sentence_generator.get_module_sentences(
        config=conf_parser, aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True,
        keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=expression_module_sentences,
                                                               module=Module.EXPRESSION)
    gene_desc.set_initial_stats(module=Module.EXPRESSION, sentence_generator=expr_sentence_generator)
    # Information poor genes
    if not gene_desc.description:
        ec_gene_id = gene_desc.gene_id[3:]
        ec_anatomy_studies = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                               expression_cluster_type=ExpressionClusterType.ANATOMY,
                                                               feature=ExpressionClusterFeature.STUDIES)
        ec_anatomy_terms = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                             feature=ExpressionClusterFeature.TERMS,
                                                             expression_cluster_type=ExpressionClusterType.ANATOMY)
        if dm.expression_ontology is not None:
            expression_enriched_module_sentences = expr_sentence_generator.get_module_sentences(
                config=conf_parser, aspect='A', qualifier="Enriched", merge_groups_with_same_prefix=True,
                keep_only_best_group=False)
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
        ec_molreg_terms = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                            expression_cluster_type=ExpressionClusterType.MOLREG,
                                                            feature=ExpressionClusterFeature.TERMS)
        ec_molereg_studies = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                               feature=ExpressionClusterFeature.STUDIES,
                                                               expression_cluster_type=ExpressionClusterType.MOLREG)

        if dm.expression_ontology is None and ec_molreg_terms:
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.EXPRESSION_CLUSTER_MOLECULE,
                description="is affected by " + concatenate_words_with_oxford_comma(ec_molreg_terms) + " based on ",
                additional_postfix_terms_list=ec_molereg_studies,
                additional_postfix_final_word="studies", use_single_form=True)


def set_do_sentence(df: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription, gene: Gene):
    do_sentence_exp_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.DO,
                                                          module=Module.DO_EXP_AND_BIO, data_manager=df,
                                                          config=conf_parser, limit_to_group="EXPERIMENTAL")
    disease_exp_module_sentences = do_sentence_exp_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_EXPERIMENTAL,
                                                               module_sentences=disease_exp_module_sentences)
    do_sentence_bio_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.DO,
                                                          module=Module.DO_EXP_AND_BIO, data_manager=df,
                                                          config=conf_parser, limit_to_group="BIOMARKER")
    disease_bio_module_sentences = do_sentence_bio_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_BIOMARKER,
                                                               module_sentences=disease_bio_module_sentences)
    do_via_orth_sentence_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.DO,
                                                               module=Module.DO_ORTHOLOGY, data_manager=df,
                                                               config=conf_parser)
    disease_via_orth_module_sentences = do_via_orth_sentence_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_ORTHOLOGY,
                                                               module_sentences=disease_via_orth_module_sentences)
    gene_desc.set_initial_stats(module=Module.DO_EXPERIMENTAL, sentence_generator=do_sentence_exp_generator)
    gene_desc.set_initial_stats(module=Module.DO_BIOMARKER, sentence_generator=do_sentence_bio_generator)
    gene_desc.set_initial_stats(module=Module.DO_ORTHOLOGY, sentence_generator=do_via_orth_sentence_generator)


def set_information_poor_sentence(orth_fullnames: List[str], selected_orthologs, ensembl_hgnc_ids_map,
                                  conf_parser: GenedescConfigParser, human_df_agr: DataManager,
                                  gene_desc: GeneDescription, dm: WBDataManager, gene: Gene):
    human_func_sent = None
    if len(orth_fullnames) == 1 and orth_fullnames[0] == "Homo sapiens":
        best_orth = get_best_human_ortholog_for_info_poor(selected_orthologs, ensembl_hgnc_ids_map,
                                                          conf_parser.get_annotations_priority(module=Module.GO),
                                                          human_df_agr, config=conf_parser)
        if best_orth:
            best_orth = "RGD:" + best_orth
            human_go_sent_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.DO,
                                                                module=Module.DO_ORTHOLOGY, data_manager=dm,
                                                                config=conf_parser, humans=True,
                                                                limit_to_group="EXPERIMENTAL")
            human_func_module_sentences = human_go_sent_generator.get_module_sentences(
                config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
            human_func_sent = human_func_module_sentences.get_description()
            if human_func_sent:
                gene_desc.set_or_extend_module_description_and_final_stats(
                    module=Module.INFO_POOR_HUMAN_FUNCTION, description="human " +
                                                                        human_df_agr.go_associations.subject_label_map[
                                                                            best_orth] + " " + human_func_sent)
    if not human_func_sent:
        protein_domains = dm.protein_domains[gene_desc.gene_id[3:]]
        if protein_domains:
            dom_word = "domain"
            if len(protein_domains) > 1:
                dom_word = "domains"
            gene_desc.set_or_extend_module_description_and_final_stats(
                module=Module.INFO_POOR_HUMAN_FUNCTION,
                description="is predicted to encode a protein with the following " + dom_word + ": " +
                            concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for
                                                                 ptdom in protein_domains]))


def set_sister_species_sentence(dm: WBDataManager, conf_parser: GenedescConfigParser, sister_sp_fullname,
                                sister_df: WBDataManager, species, organism, gene_desc: GeneDescription, gene: Gene):
    best_ortholog = dm.get_best_orthologs_for_gene(
        gene_desc.gene_id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
        ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                             "HEP"])[0][0]
    sister_sentences_generator = OntologySentenceGenerator(gene_id=gene.id, annot_type=DataType.DO,
                                                           module=Module.DO_ORTHOLOGY, data_manager=dm,
                                                           config=conf_parser,
                                                           humans=sister_sp_fullname == "Homo sapiens",
                                                           limit_to_group="EXPERIMENTAL")
    sister_sp_module_sentences = sister_sentences_generator.get_module_sentences(
        config=conf_parser, aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True).get_description()
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
    for organism in organisms_list:
        logger.info("Processing organism " + organism)
        species = conf_parser.get_wb_organisms_info()
        dm, sister_df, df_agr, orth_fullnames, sister_sp_fullname = load_data(
            species=species, organism=organism, conf_parser=conf_parser)
        desc_writer = DescriptionsWriter()
        desc_writer.overall_properties.species = organism
        desc_writer.overall_properties.release_version = conf_parser.get_wb_release()[0:-1] + str(
            int(conf_parser.get_wb_release()[-1]) + 1)
        desc_writer.overall_properties.date = datetime.date.today().strftime("%B %d, %Y")
        for gene in dm.get_gene_data():
            logger.debug("Generating description for gene " + gene.name)

            gene_desc = GeneDescription(gene_id=gene.id, gene_name=gene.name, add_gene_name=True)
            selected_orthologs = set_orthology_sentence(dm=dm, orth_fullnames=orth_fullnames,
                                                        human_genes_props=human_genes_props, gene_desc=gene_desc,
                                                        tpc_token=args.textpresso_token)
            set_go_sentences(dm=dm, conf_parser=conf_parser, gene_desc=gene_desc, gene=gene)
            set_expression_sentence(dm=dm, conf_parser=conf_parser, gene_desc=gene_desc, gene=gene)
            set_do_sentence(df=dm, conf_parser=conf_parser, gene=gene, gene_desc=gene_desc)
            if not gene_desc.go_description:
                set_information_poor_sentence(orth_fullnames=orth_fullnames,
                                              selected_orthologs=selected_orthologs,
                                              ensembl_hgnc_ids_map=ensembl_hgnc_ids_map, conf_parser=conf_parser,
                                              human_df_agr=df_agr, gene_desc=gene_desc, dm=dm, gene=gene)
            if "main_sister_species" in species[organism] and species[organism]["main_sister_species"] and \
                    dm.get_best_orthologs_for_gene(gene.id, orth_species_full_name=[sister_sp_fullname],
                                                   sister_species_data_fetcher=sister_df,
                                                   ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP",
                                                                        "HDA", "HMP", "HGI", "HEP"])[0]:
                set_sister_species_sentence(dm=dm, sister_sp_fullname=sister_sp_fullname, sister_df=sister_df,
                                            species=species, organism=organism, gene_desc=gene_desc,
                                            conf_parser=conf_parser, gene=gene)
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
