#!/usr/bin/env python3

import argparse
import datetime

from genedescriptions.data_fetcher import WBDataFetcher, DataFetcher
from genedescriptions.descriptions_rules import *
from genedescriptions.descriptions_writer import JsonGDWriter


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
    parser.add_argument("-v", "--output-version", metavar="version_number", dest="version_number", type=str,
                        help="release version number")

    args = parser.parse_args()
    conf_parser = GenedescConfigParser(args.config_file)
    logging.basicConfig(filename=args.log_file, level=args.log_level, format='%(asctime)s - %(name)s - %(levelname)s:'
                                                                             '%(message)s')

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
                            "add_multiple_if_covers_more_children": False}
    do_sent_gen_common_props = {"evidence_groups_priority_list": conf_parser.get_do_evidence_groups_priority_list(),
                                "prepostfix_sentences_map": conf_parser.get_do_prepostfix_sentences_map(),
                                "prepostfix_special_cases_sent_map": None,
                                "evidence_codes_groups_map": conf_parser.get_do_evidence_codes_groups_map()}
    do_sent_common_props = {"remove_parent_terms": conf_parser.get_do_remove_parents_if_children_are_present(),
                            "remove_child_terms": conf_parser.get_do_remove_children_if_parent_is_present(),
                            "merge_num_terms_threshold": conf_parser.get_do_trim_min_num_terms(),
                            "merge_max_num_terms": conf_parser.get_do_max_num_terms(),
                            "merge_min_distance_from_root": conf_parser.get_do_trim_min_distance_from_root(),
                            "truncate_others_generic_word": conf_parser.get_do_truncate_others_aggregation_word(),
                            "truncate_others_aspect_words": conf_parser.get_do_truncate_others_terms(),
                            "add_multiple_if_covers_more_children": True}
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
        "add_multiple_if_covers_more_children": False, "rename_cell": True}

    organisms_list = conf_parser.get_wb_organisms_to_process()
    human_genes_props = DataFetcher.get_human_gene_props()
    ensembl_hgnc_ids_map = DataFetcher.get_ensembl_hgnc_ids_map()
    for organism in organisms_list:
        logging.info("Processing organism " + organism)
        sister_df = None
        species = conf_parser.get_wb_species()
        sister_sp_fullname = ""
        df_agr = None
        if "main_sister_species" in species[organism] and "full_name" in \
                species[species[organism]["main_sister_species"]]:
            sister_sp_fullname = species[species[organism]["main_sister_species"]]["full_name"]
        orthologs_sp_fullname = ""
        if "ortholog" in species[organism] and all(["full_name" in species[ortholog_sp] for ortholog_sp in
                                                    species[organism]["ortholog"]]):
            orthologs_sp_fullname = [species[ortholog_sp]["full_name"] for ortholog_sp in species[organism]["ortholog"]]
        df = WBDataFetcher(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                           release_version=conf_parser.get_release("wb_data_fetcher"),
                           species=organism, project_id=species[organism]["project_id"],
                           cache_location=conf_parser.get_cache_location(), do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"], sister_sp_fullname=sister_sp_fullname)
        if organism == "c_elegans":
            df_agr = DataFetcher(go_relations=["subClassOf", "BFO:0000050"], do_relations=None)
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
            sister_df = WBDataFetcher(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                                      release_version=conf_parser.get_release("wb_data_fetcher"),
                                      species=species[organism]["main_sister_species"],
                                      project_id=species[species[organism]["main_sister_species"]]["project_id"],
                                      cache_location=conf_parser.get_cache_location(), do_relations=None,
                                      go_relations=["subClassOf", "BFO:0000050"])
            logging.info("Loading all data for sister species")
            sister_df.load_all_data_from_file(go_terms_replacement_regex=conf_parser.get_go_rename_terms(),
                                              go_terms_exclusion_list=conf_parser.get_go_terms_exclusion_list(),
                                              do_terms_replacement_regex=None,
                                              do_terms_exclusion_list=conf_parser.get_do_terms_exclusion_list())
        logging.info("Loading all data for main species")
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
        desc_writer = JsonGDWriter()
        desc_writer.overall_properties.species = organism
        desc_writer.overall_properties.release_version = conf_parser.get_release("wb_data_fetcher")
        desc_writer.overall_properties.date = datetime.date.today().strftime("%B %d, %Y")
        for gene in df.get_gene_data():
            logging.debug("Generating description for gene " + gene.name)

            compose_wormbase_description(gene=gene, conf_parser=conf_parser, species=species, organism=organism, df=df,
                                         orthologs_sp_fullname=orthologs_sp_fullname,
                                         go_sent_gen_common_props=go_sent_gen_common_props,
                                         go_sent_common_props=go_sent_common_props, human_genes_props=human_genes_props,
                                         do_sent_gen_common_prop=do_sent_gen_common_props,
                                         do_sent_common_props=do_sent_common_props,
                                         sister_sp_fullname=sister_sp_fullname, sister_df=sister_df,
                                         human_df_agr=df_agr, desc_writer=desc_writer,
                                         ensembl_hgnc_ids_map=ensembl_hgnc_ids_map,
                                         expr_sent_gen_common_props=expr_sent_gen_common_props,
                                         expr_sent_common_props=expr_sent_common_props)
        desc_writer.write(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                       organism + "_with_stats.json"), pretty=True, include_single_gene_stats=True)
        desc_writer.write(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                       organism + "_no_stats.json"), pretty=True, include_single_gene_stats=False)


if __name__ == '__main__':
    main()
