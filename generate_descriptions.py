#!/usr/bin/env python3

import argparse
import configparser
from data_fetcher import WBRawDataFetcher
from descriptions_rules import *
import logging


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="genedesc.ini", help="configuration file")
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

    config = configparser.ConfigParser()
    config.read(args.config_file)

    logging.basicConfig(filename=args.log_file, level=args.log_level)

    raw_files_source = config.get("data_fetcher", "raw_files_source")
    chebi_files_source = config.get("data_fetcher", "chebi_files_source")
    cache_location = config.get("data_fetcher", "cache_location")
    species = config.get("generic", "species").split(",")
    project_ids = config.get("generic", "project_ids").split(",")
    evidence_codes = config.get("go_sentences_options", "evidence_codes").split(",")
    evidence_codes = dict(zip(evidence_codes, range(len(evidence_codes))))
    evidence_groups = config.get("go_sentences_options", "evidence_groups").split(",")
    evidence_groups = dict(zip(evidence_groups, range(len(evidence_groups))))
    evidence_codes_groups_map = config.get("go_sentences_options", "evidence_codes_groups_map").split(",")
    evidence_codes_groups_map = dict(zip(evidence_codes, [int(num) for num in evidence_codes_groups_map]))
    go_prepostfix_sentences_list = config.get("go_sentences_options", "go_prepostfix_sentences_map").split(";")
    go_prepostfix_sentences_map = {}
    for go_prepostfix_sentence in go_prepostfix_sentences_list:
        indices, values = go_prepostfix_sentence.split(":")
        aspect, group = indices.split(",")
        prefix, postfix = values.split(",")
        go_prepostfix_sentences_map[(GO_ASPECT[aspect], evidence_groups[group])] = (prefix, postfix)
    go_annotations_priority = config.get("go_sentences_options", "go_annotations_priority").split(",")

    df = WBRawDataFetcher(raw_files_source=raw_files_source, chebi_file_url=chebi_files_source,
                          release_version=args.wormbase_number, species=species[3],
                          project_id=project_ids[3], cache_location=cache_location, use_cache=args.use_cache)

    df.load_go_data()
    for gene in df.get_gene_data():
        print(gene.id, gene.name)
        sentences = generate_go_sentences(df.get_go_annotations(gene.id, priority_list=go_annotations_priority),
                                          evidence_groups, go_prepostfix_sentences_map,
                                          evidence_codes_groups_map)
        if sentences:
            joined_sent = []
            func_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                GO_ASPECT["MOLECULAR_FUNCTION"])])
            if func_sent:
                joined_sent.append(func_sent)
            proc_sent = " and ".join([sentence.text for sentence in sentences.get_sentences(
                GO_ASPECT["BIOLOGICAL_PROCESS"])])
            if proc_sent:
                joined_sent.append(proc_sent)
            comp_sent = "and ".join([sentence.text for sentence in sentences.get_sentences(
                GO_ASPECT["CELLULAR_COMPONENT"])])
            if comp_sent:
                joined_sent.append(comp_sent)

            go_desc = "; ".join(joined_sent) + "."
            print(go_desc.capitalize())
        else:
            print("No description available")
        print()


if __name__ == '__main__':
    main()
