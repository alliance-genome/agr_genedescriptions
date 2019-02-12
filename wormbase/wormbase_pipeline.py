#!/usr/bin/env python3

import argparse
import datetime
import logging
import os

from typing import List
from num2words import num2words

from genedescriptions.api_manager import APIManager
from genedescriptions.commons import DataType, Module, Gene
from genedescriptions.config_parser import GenedescConfigParser, ConfigModuleProperty
from genedescriptions.data_manager import DataManager, ExpressionClusterType, ExpressionClusterFeature
from genedescriptions.gene_description import GeneDescription
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.descriptions_writer import DescriptionsWriter
from genedescriptions.precanned_modules import set_gene_ontology_module, set_disease_module, \
    generate_ortholog_sentence_wormbase_human, generate_ortholog_sentence_wormbase_non_c_elegans
from genedescriptions.sentence_generation_functions import concatenate_words_with_oxford_comma, \
    get_best_human_ortholog_for_info_poor
from wormbase.wb_data_manager import WBDataManager


def load_data(organism, conf_parser: GenedescConfigParser):
    logger = logging.getLogger("WB Gene Description Pipeline - Data loader")
    sister_df = None
    df_agr = None
    organisms_info = conf_parser.get_wb_organisms_info()
    df = WBDataManager(species=organism, do_relations=None, go_relations=["subClassOf", "BFO:0000050"],
                       config=conf_parser)
    if organism == "c_elegans":
        df_agr = DataManager(go_relations=["subClassOf", "BFO:0000050"], do_relations=None)
        df_agr.load_ontology_from_file(ontology_type=DataType.GO,
                                       ontology_url=conf_parser.get_wb_human_orthologs_go_ontology(),
                                       ontology_cache_path=os.path.join(conf_parser.get_cache_dir(),
                                                                        "wormbase_agr_human", "go_ontology.obo"),
                                       config=conf_parser)
        df_agr.load_associations_from_file(associations_type=DataType.GO,
                                           associations_url=conf_parser.get_wb_human_orthologs_go_associations(),
                                           associations_cache_path=os.path.join(
                                               conf_parser.get_cache_dir(), "wormbase_agr_human", "go_assoc.daf.gz"),
                                           config=conf_parser)
    if "main_sister_species" in organisms_info[organism] and organisms_info[organism]["main_sister_species"]:
        sister_df = WBDataManager(species=organisms_info[organism]["main_sister_species"],
                                  do_relations=None, go_relations=["subClassOf", "BFO:0000050"], config=conf_parser)
        logger.info("Loading GO data for sister species")
        sister_df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url=sister_df.go_ontology_url,
                                          ontology_cache_path=sister_df.go_ontology_cache_path,
                                          config=conf_parser)
        sister_df.load_associations_from_file(associations_type=DataType.GO,
                                              associations_url=sister_df.go_associations_url,
                                              associations_cache_path=sister_df.go_associations_cache_path,
                                              config=conf_parser)
    logger.info("Loading all data for main species")
    df.load_all_data_from_file()
    return df, sister_df, df_agr


def set_orthology_sentence(dm: WBDataManager, orth_fullnames: List[str], gene_desc: GeneDescription,
                           human_genes_props, api_manager):
    best_orthologs, selected_orth_name = dm.get_best_orthologs_for_gene(gene_desc.gene_id,
                                                                        orth_species_full_name=orth_fullnames)
    selected_orthologs = []
    if best_orthologs:
        gene_desc.stats.set_best_orthologs = [orth[0] for orth in best_orthologs]
        if len(orth_fullnames) == 1 and orth_fullnames[0] == "Homo sapiens":
            sel_orthologs, orth_sent = generate_ortholog_sentence_wormbase_human(best_orthologs, human_genes_props)
            selected_orthologs = [orth for orth in best_orthologs if orth[1] in sel_orthologs]
        else:
            orth_sent = generate_ortholog_sentence_wormbase_non_c_elegans(best_orthologs, selected_orth_name,
                                                                          api_manager=api_manager)
        gene_desc.set_or_extend_module_description_and_final_stats(module=Module.ORTHOLOGY, description=orth_sent)
    return selected_orthologs


def set_tissue_expression_sentence(dm, gene, conf_parser, gene_desc):
    expr_sentence_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.EXPRESSION, data_manager=dm,
                                                        config=conf_parser)
    expression_module_sentences = expr_sentence_generator.get_module_sentences(
        config=conf_parser, aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True,
        keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=expression_module_sentences,
                                                               module=Module.EXPRESSION)
    gene_desc.set_initial_stats(module=Module.EXPRESSION, sentence_generator=expr_sentence_generator)


def set_expression_cluster_sentence(dm: WBDataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription,
                                    gene: Gene, api_manager: APIManager):

    expr_sentence_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.EXPRESSION, data_manager=dm,
                                                        config=conf_parser)
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
            module=Module.EXPRESSION_CLUSTER_ANATOMY,
            description=expression_enriched_module_sentences.get_description(),
            additional_postfix_terms_list=ec_anatomy_studies, additional_postfix_final_word="studies",
            use_single_form=True)
    elif ec_anatomy_terms:
        gene_desc.set_or_extend_module_description_and_final_stats(
            module=Module.EXPRESSION_CLUSTER_ANATOMY,
            description="is enriched in " + concatenate_words_with_oxford_comma(ec_anatomy_terms) + " based on",
            additional_postfix_terms_list=ec_anatomy_studies,
            additional_postfix_final_word="studies", use_single_form=True)
    ec_molreg_terms = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                        expression_cluster_type=ExpressionClusterType.MOLREG,
                                                        feature=ExpressionClusterFeature.TERMS)
    ec_molreg_studies = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                          feature=ExpressionClusterFeature.STUDIES,
                                                          expression_cluster_type=ExpressionClusterType.MOLREG)
    ec_genereg_terms = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                         expression_cluster_type=ExpressionClusterType.GENEREG,
                                                         feature=ExpressionClusterFeature.TERMS)
    ec_genereg_studies = dm.get_expression_cluster_feature(gene_id=ec_gene_id,
                                                           feature=ExpressionClusterFeature.STUDIES,
                                                           expression_cluster_type=ExpressionClusterType.GENEREG)
    if ec_genereg_terms:
        several_word = ""
        if len(ec_genereg_terms) > 3:
            t_p = [t_p for t_p in sorted([[term, api_manager.get_textpresso_popularity(term)] for
                                          term in ec_genereg_terms], key=lambda x: (x[1], x[0][1]),
                                         reverse=True)]
            ec_genereg_terms = [term for term, popularity in t_p[0:3]]
            several_word = "several genes including "
        gene_desc.set_or_extend_module_description_and_final_stats(
            module=Module.EXPRESSION_CLUSTER_GENE,
            description="is affected by " + several_word +
                        concatenate_words_with_oxford_comma(ec_genereg_terms) + " based on",
            additional_postfix_terms_list=ec_genereg_studies,
            additional_postfix_final_word="studies", use_single_form=True)
    if ec_molreg_terms:
        several_word = ""
        if len(ec_molreg_terms) > 3:
            several_word = num2words(len(ec_molreg_terms)) + " chemicals including "
        gene_desc.set_or_extend_module_description_and_final_stats(
            module=Module.EXPRESSION_CLUSTER_MOLECULE,
            description="is affected by " + several_word + concatenate_words_with_oxford_comma(
                ec_molreg_terms[0:3]) + " based on",
            additional_postfix_terms_list=ec_molreg_studies,
            additional_postfix_final_word="studies", use_single_form=True)


def set_information_poor_sentence(orth_fullnames: List[str], selected_orthologs, ensembl_hgnc_ids_map,
                                  conf_parser: GenedescConfigParser, human_df_agr: DataManager,
                                  gene_desc: GeneDescription, dm: WBDataManager, gene: Gene):
    if len(orth_fullnames) == 1 and orth_fullnames[0] == "Homo sapiens":
        best_orth = get_best_human_ortholog_for_info_poor(selected_orthologs, ensembl_hgnc_ids_map,
                                                          conf_parser.get_annotations_priority(module=Module.GO),
                                                          human_df_agr, config=conf_parser)
        if best_orth:
            if not best_orth.startswith("RGD:"):
                best_orth = "RGD:" + best_orth
            human_go_sent_generator = OntologySentenceGenerator(gene_id=best_orth, module=Module.GO,
                                                                data_manager=human_df_agr, config=conf_parser,
                                                                humans=False, limit_to_group="EXPERIMENTAL")
            human_func_module_sentences = human_go_sent_generator.get_module_sentences(
                config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
            human_func_sent = human_func_module_sentences.get_description()
            if human_func_sent:
                gene_desc.set_or_extend_module_description_and_final_stats(
                    module=Module.INFO_POOR_HUMAN_FUNCTION, description="human " +
                                                                        human_df_agr.go_associations.subject_label_map[
                                                                            best_orth] + " " + human_func_sent)

    protein_domains = dm.protein_domains[gene_desc.gene_id[3:]]
    if protein_domains:
        dom_word = "domain"
        if len(protein_domains) > 1:
            dom_word = "domains"
        gene_desc.set_or_extend_module_description_and_final_stats(
            module=Module.PROTEIN_DOMAIN,
            description="is predicted to encode a protein with the following " + dom_word + ": " +
                        concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for
                                                             ptdom in protein_domains]))


def set_sister_species_sentence(dm: WBDataManager, conf_parser: GenedescConfigParser, sister_sp_fullname,
                                sister_df: WBDataManager, species, organism, gene_desc: GeneDescription, gene: Gene):
    best_ortholog = dm.get_best_orthologs_for_gene(
        gene_desc.gene_id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
        ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                             "HEP"])[0][0]
    if not best_ortholog[0].startswith("WB:"):
        best_ortholog[0] = "WB:" + best_ortholog[0]
    sister_sentences_generator = OntologySentenceGenerator(gene_id=best_ortholog[0], module=Module.GO,
                                                           data_manager=sister_df, config=conf_parser,
                                                           humans=sister_sp_fullname == "Homo sapiens",
                                                           limit_to_group="EXPERIMENTAL")
    sister_sp_module_sentences = sister_sentences_generator.get_module_sentences(
        config=conf_parser, aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True)
    if sister_sp_module_sentences.contains_sentences():
        gene_desc.set_or_extend_module_description_and_final_stats(
            module=Module.SISTER_SP, description="in " + species[species[organism]["main_sister_species"]]["name"] +
                                                 ", " + best_ortholog[1] + " " +
                                                 sister_sp_module_sentences.get_description())


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
    api_manager = APIManager(textpresso_api_token=args.textpresso_token)
    for organism in organisms_list:
        logger.info("Processing organism " + organism)
        species = conf_parser.get_wb_organisms_info()
        dm, sister_df, df_agr = load_data(organism=organism, conf_parser=conf_parser)
        desc_writer = DescriptionsWriter()
        desc_writer.overall_properties.species = organism
        desc_writer.overall_properties.release_version = conf_parser.get_wb_release()[0:-1] + str(
            int(conf_parser.get_wb_release()[-1]) + 1)
        desc_writer.overall_properties.date = datetime.date.today().strftime("%B %d, %Y")
        for gene in dm.get_gene_data():
            logger.debug("Generating description for gene " + gene.name)
            gene_desc = GeneDescription(gene_id=gene.id, gene_name=gene.name, add_gene_name=False)
            selected_orthologs = set_orthology_sentence(dm=dm, orth_fullnames=dm.orth_fullnames,
                                                        human_genes_props=human_genes_props, gene_desc=gene_desc,
                                                        api_manager=api_manager)
            set_gene_ontology_module(dm=dm, conf_parser=conf_parser, gene_desc=gene_desc, gene=gene)
            set_tissue_expression_sentence(dm=dm, gene=gene, conf_parser=conf_parser, gene_desc=gene_desc)
            if not gene_desc.description:
                set_expression_cluster_sentence(dm=dm, conf_parser=conf_parser, gene_desc=gene_desc, gene=gene,
                                                api_manager=api_manager)
            set_disease_module(df=dm, conf_parser=conf_parser, gene=gene, gene_desc=gene_desc)
            if not gene_desc.go_description:
                set_information_poor_sentence(orth_fullnames=dm.orth_fullnames,
                                              selected_orthologs=selected_orthologs,
                                              ensembl_hgnc_ids_map=ensembl_hgnc_ids_map, conf_parser=conf_parser,
                                              human_df_agr=df_agr, gene_desc=gene_desc, dm=dm, gene=gene)
            if "main_sister_species" in species[organism] and species[organism]["main_sister_species"] and \
                    dm.get_best_orthologs_for_gene(gene.id, orth_species_full_name=[dm.sister_sp_fullname],
                                                   sister_species_data_fetcher=sister_df,
                                                   ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP",
                                                                        "HDA", "HMP", "HGI", "HEP"])[0]:
                set_sister_species_sentence(dm=dm, sister_sp_fullname=dm.sister_sp_fullname, sister_df=sister_df,
                                            species=species, organism=organism, gene_desc=gene_desc,
                                            conf_parser=conf_parser, gene=gene)
            desc_writer.add_gene_desc(gene_desc)
        logger.info("All genes processed for " + organism)
        date_prefix = datetime.date.today().strftime("%Y%m%d")
        if "json" in args.output_formats:
            logger.info("Writing descriptions to json")
            desc_writer.write_json(os.path.join(conf_parser.get_out_dir(), date_prefix + "_" + organism + ".json"),
                                   pretty=True, include_single_gene_stats=True, data_manager=dm)
        if "txt" in args.output_formats:
            logger.info("Writing descriptions to txt")
            desc_writer.write_plain_text(os.path.join(conf_parser.get_out_dir(), date_prefix + "_" + organism + ".txt"))
        if "tsv" in args.output_formats:
            logger.info("Writing descriptions to tsv")
            desc_writer.write_tsv(os.path.join(conf_parser.get_out_dir(), date_prefix + "_" + organism + ".tsv"))
        if "ace" in args.output_formats:
            logger.info("Writing descriptions to ace")
            curators = ["WBPerson324", "WBPerson37462"]
            release_version = conf_parser.get_wb_release()
            desc_writer.write_ace(os.path.join(conf_parser.get_out_dir(), date_prefix + "_" + organism + ".ace"),
                                  curators, release_version)


if __name__ == '__main__':
    main()
