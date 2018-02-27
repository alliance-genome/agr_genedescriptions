#!/usr/bin/env python3

import argparse
import json
import logging

from config_parser import GenedescConfigParser
from data_fetcher import WBRawDataFetcher, AGRRawDataFetcher
from descriptions_rules import *


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.json", help="configuration file")
    parser.add_argument("-C", "--use-cache", dest="use_cache", action="store_true", default=False,
                        help="Use cached source files from cache_location specified in config file. Download them from "
                             "raw_file_source (configured in config file) if not yet cached")
    parser.add_argument("-l", "--log-file", metavar="log_file", dest="log_file", type=str,
                        default="genedescriptions.log", help="path to the log file to generate")
    parser.add_argument("-L", "--log-level", dest="log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], help="set the logging level")
    parser.add_argument("-v", "--output-version", metavar="version_number", dest="version_number", type=str,
                        help="release version number")

    args = parser.parse_args()

    conf_parser = GenedescConfigParser(args.config_file)

    logging.basicConfig(filename=args.log_file, level=args.log_level)

    cache_location = conf_parser.get_cache_location()
    species = conf_parser.get_wb_species()
    go_prepostfix_sentences_map = conf_parser.get_go_prepostfix_sentences_map()
    go_prepostfix_special_cases_sent_map = conf_parser.get_go_prepostfix_special_cases_sent_map()
    go_annotations_priority = conf_parser.get_go_annotations_priority()
    evidence_groups_priority_list = conf_parser.get_evidence_groups_priority_list()
    evidence_codes_groups_map = conf_parser.get_evidence_codes_groups_map()
    go_terms_exclusion_list = conf_parser.get_go_terms_exclusion_list()

    if conf_parser.get_data_fetcher() == "agr_data_fetcher":
        organisms_list = conf_parser.get_agr_organisms_to_process()
    else:
        organisms_list = conf_parser.get_wb_organisms_to_process()
    df = None
    for organism in organisms_list:
        if conf_parser.get_data_fetcher() == "agr_data_fetcher":
            df = AGRRawDataFetcher(go_terms_exclusion_list=go_terms_exclusion_list,
                                   go_terms_replacement_dict=conf_parser.get_go_rename_terms(),
                                   raw_files_source=conf_parser.get_raw_file_sources(conf_parser.get_data_fetcher()),
                                   chebi_file_url=conf_parser.get_chebi_file_source(),
                                   release_version=conf_parser.get_release(conf_parser.get_data_fetcher()),
                                   main_file_name=conf_parser.get_agr_mod_property(organism, "main_files"),
                                   bgi_file_name=conf_parser.get_agr_mod_property(organism, "bgi_file"),
                                   go_annotations_file_name=conf_parser.get_agr_mod_property(organism,
                                                                                             "go_annotations"),
                                   organism_name=conf_parser.get_agr_mod_property(organism, "name"),
                                   cache_location=cache_location, use_cache=args.use_cache)
        else:
            df = WBRawDataFetcher(go_terms_exclusion_list=go_terms_exclusion_list,
                                  go_terms_replacement_dict=conf_parser.get_go_rename_terms(),
                                  raw_files_source=conf_parser.get_raw_file_sources(conf_parser.get_data_fetcher()),
                                  chebi_file_url=conf_parser.get_chebi_file_source(),
                                  release_version=conf_parser.get_release(conf_parser.get_data_fetcher()),
                                  species=organism,
                                  project_id=species[organism]["project_id"],
                                  cache_location=cache_location, use_cache=args.use_cache)

    df.load_go_data()
    json_gene_desc = []
    for gene in df.get_gene_data():
        gene_dict = {"name": gene.name}
        sentences = generate_go_sentences(df.get_go_annotations(gene.id, priority_list=go_annotations_priority),
                                          go_ontology=df.get_go_ontology(),
                                          evidence_groups_priority_list=evidence_groups_priority_list,
                                          go_prepostfix_sentences_map=go_prepostfix_sentences_map,
                                          go_prepostfix_special_cases_sent_map=go_prepostfix_special_cases_sent_map,
                                          evidence_codes_groups_map=evidence_codes_groups_map,
                                          remove_parent_terms=True, merge_num_terms_threshold=3,
                                          merge_min_distance_from_root=2)
        num_terms_f = 0
        num_terms_p = 0
        num_terms_c = 0
        if sentences:
            joined_sent = []
            func_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='F', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True)])
            if func_sent:
                joined_sent.append(func_sent)
                num_terms_f = sum(map(len, [sentence.terms for sentence in sentences.get_sentences(
                    go_aspect='F', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True)]))
            proc_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='P', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True,
                keep_only_best_group=True)])
            if proc_sent:
                joined_sent.append(proc_sent)
                num_terms_p = sum(map(len, [sentence.terms for sentence in sentences.get_sentences(
                    go_aspect='P', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True,
                    keep_only_best_group=True)]))
            comp_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='C', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True,
                keep_only_best_group=True)])
            if comp_sent:
                joined_sent.append(comp_sent)
                num_terms_c = sum(map(len, [sentence.terms for sentence in sentences.get_sentences(
                    go_aspect='C', go_ontology=df.get_go_ontology(), merge_groups_with_same_prefix=True,
                    keep_only_best_group=True)]))

            go_desc = "; ".join(joined_sent) + "."
            gene_dict["description"] = go_desc.capitalize()
        else:
            gene_dict["description"] = "No description available"
        gene_dict["stats"] = {"num_terms_f": num_terms_f, "num_terms_p": num_terms_p, "num_terms_c": num_terms_c }
        json_gene_desc.append(gene_dict)

    json_string = json.dumps(json_gene_desc)
    print(json_string)


if __name__ == '__main__':
    main()
