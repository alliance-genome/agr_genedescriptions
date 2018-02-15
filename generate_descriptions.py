#!/usr/bin/env python3

import argparse
import json

from config_parser import GenedescConfigParser
from data_fetcher import WBRawDataFetcher, AGRRawDataFetcher
from descriptions_rules import *
import logging


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

    if conf_parser.get_data_fetcher() == "AGR":
        df = AGRRawDataFetcher(go_terms_exclusion_list=go_terms_exclusion_list,
                               raw_files_source=conf_parser.get_raw_file_sources(conf_parser.get_data_fetcher()),
                               chebi_file_url=conf_parser.get_chebi_file_source(),
                               release_version=conf_parser.get_release(conf_parser.get_data_fetcher()),
                               gff_file_name=conf_parser.get_agr_mod_property("zfin", "gff"),
                               go_annotations_file_name=conf_parser.get_agr_mod_property("zfin", "go_annotations"),
                               organism_name=conf_parser.get_agr_mod_property("zfin", "name"),
                               cache_location=cache_location, use_cache=args.use_cache)
    else:
        df = WBRawDataFetcher(go_terms_exclusion_list=go_terms_exclusion_list,
                              raw_files_source=conf_parser.get_raw_file_sources(conf_parser.get_data_fetcher()),
                              chebi_file_url=conf_parser.get_chebi_file_source(),
                              release_version=conf_parser.get_release(conf_parser.get_data_fetcher()),
                              species="c_elegans",
                              project_id=species["c_elegans"]["project_id"],
                              cache_location=cache_location, use_cache=args.use_cache)

    df.load_go_data()
    for gene in df.get_gene_data():
        print(gene.id, gene.name)
        sentences = generate_go_sentences(df.get_go_annotations(gene.id, priority_list=go_annotations_priority),
                                          evidence_groups_priority_list=evidence_groups_priority_list,
                                          go_prepostfix_sentences_map=go_prepostfix_sentences_map,
                                          go_prepostfix_special_cases_sent_map=go_prepostfix_special_cases_sent_map,
                                          evidence_codes_groups_map=evidence_codes_groups_map)
        if sentences:
            joined_sent = []
            func_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='F', merge_groups_with_same_prefix=True)])
            if func_sent:
                joined_sent.append(func_sent)
            proc_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True)])
            if proc_sent:
                joined_sent.append(proc_sent)
            comp_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                go_aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True)])
            if comp_sent:
                joined_sent.append(comp_sent)

            go_desc = "; ".join(joined_sent) + "."
            print(go_desc.capitalize())
        else:
            print("No description available")
        print()


if __name__ == '__main__':
    main()
