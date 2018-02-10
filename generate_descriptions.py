#!/usr/bin/env python3

import argparse
import json

from data_fetcher import WBRawDataFetcher
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
    parser.add_argument("-w", "--wormbase-version", metavar="wormbase_number", dest="wormbase_number", type=str,
                        help="wormbase input files version number")

    args = parser.parse_args()

    with open(args.config_file) as conf_file:
        config = json.load(conf_file)

    logging.basicConfig(filename=args.log_file, level=args.log_level)

    raw_files_source = config["wb_data_fetcher_options"]["raw_files_source"]
    chebi_file_source = config["wb_data_fetcher_options"]["chebi_file_source"]
    cache_location = config["generic_data_fetcher_options"]["cache_location"]
    species = config["wb_data_fetcher_options"]["species"]
    evidence_codes = config["go_sentences_options"]["evidence_codes"]
    go_prepostfix_sentences_map = {(prepost["aspect"], prepost["group"]): (prepost["prefix"], prepost["postfix"]) for
                                   prepost in config["go_sentences_options"]["go_prepostfix_sentences_map"]}
    go_annotations_priority = [name for name, priority in sorted([(ec["name"], ec["priority"]) for ec
                                                                  in evidence_codes], key=lambda x: x[1])]
    evidence_groups_list = list(set([evidence["group"] for evidence in evidence_codes]))
    evidence_codes_groups_map = {evidence["name"]: evidence["group"] for evidence in evidence_codes}

    df = WBRawDataFetcher(raw_files_source=raw_files_source, chebi_file_url=chebi_file_source,
                          release_version=args.wormbase_number, species=species[3]["name"],
                          project_id=species[3]["project_id"], cache_location=cache_location, use_cache=args.use_cache)

    df.load_go_data()
    for gene in df.get_gene_data():
        print(gene.id, gene.name)
        sentences = generate_go_sentences(df.get_go_annotations(gene.id, priority_list=go_annotations_priority),
                                          evidence_groups_list, go_prepostfix_sentences_map,
                                          evidence_codes_groups_map)
        if sentences:
            joined_sent = []
            func_sent = " and ".join([sentence.text for sentence in sentences.get_sentences('F')])
            if func_sent:
                joined_sent.append(func_sent)
            proc_sent = " and ".join([sentence.text for sentence in sentences.get_sentences('B')])
            if proc_sent:
                joined_sent.append(proc_sent)
            comp_sent = "and ".join([sentence.text for sentence in sentences.get_sentences('C')])
            if comp_sent:
                joined_sent.append(comp_sent)

            go_desc = "; ".join(joined_sent) + "."
            print(go_desc.capitalize())
        else:
            print("No description available")
        print()


if __name__ == '__main__':
    main()
