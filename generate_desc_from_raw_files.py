#!/usr/bin/env python3

import argparse
import os

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_fetcher import WBDataFetcher, DataType
from genedescriptions.descriptions_rules import *
from genedescriptions.descriptions_writer import JsonGDWriter, GeneDesc


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.yml", help="configuration file. Default ./config.yaml")
    parser.add_argument("-C", "--use-cache", dest="use_cache", action="store_true", default=False,
                        help="Use cached source files from cache_location specified in config file. Download them from "
                             "raw_file_source (configured in config file) if not yet cached")
    parser.add_argument("-l", "--log-file", metavar="log_file", dest="log_file", type=str,
                        default="genedescriptions.log",
                        help="path to the log file to generate. Default ./genedescriptions.log")
    parser.add_argument("-L", "--log-level", dest="log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], help="set the logging level")
    parser.add_argument("-v", "--output-version", metavar="version_number", dest="version_number", type=str,
                        help="release version number")

    args = parser.parse_args()
    conf_parser = GenedescConfigParser(args.config_file)
    logging.basicConfig(filename=args.log_file, level=args.log_level)

    go_sent_gen_common_props = {"evidence_groups_priority_list": conf_parser.get_go_evidence_groups_priority_list(),
                                "prepostfix_sentences_map": conf_parser.get_go_prepostfix_sentences_map(),
                                "prepostfix_special_cases_sent_map": conf_parser.get_go_prepostfix_special_cases_sent_map(),
                                "evidence_codes_groups_map": conf_parser.get_go_evidence_codes_groups_map()}
    go_sent_common_props = {"remove_parent_terms": conf_parser.get_go_remove_parents_if_children_are_present(),
                            "merge_num_terms_threshold": conf_parser.get_go_trim_min_num_terms(),
                            "merge_min_distance_from_root": conf_parser.get_go_trim_min_distance_from_root(),
                            "truncate_others_generic_word": conf_parser.get_go_truncate_others_aggregation_word(),
                            "truncate_others_aspect_words": conf_parser.get_go_truncate_others_terms(),
                            "add_multiple_if_covers_more_children": False}
    do_sent_gen_common_prop = {"evidence_groups_priority_list": conf_parser.get_do_evidence_groups_priority_list(),
                               "prepostfix_sentences_map": conf_parser.get_do_prepostfix_sentences_map(),
                               "prepostfix_special_cases_sent_map": None,
                               "evidence_codes_groups_map": conf_parser.get_do_evidence_codes_groups_map()}
    do_sent_common_props = {"remove_parent_terms": conf_parser.get_do_remove_parents_if_children_are_present(),
                            "merge_num_terms_threshold": conf_parser.get_do_trim_min_num_terms(),
                            "merge_min_distance_from_root": conf_parser.get_do_trim_min_distance_from_root(),
                            "truncate_others_generic_word": conf_parser.get_do_truncate_others_aggregation_word(),
                            "truncate_others_aspect_words": conf_parser.get_do_truncate_others_terms(),
                            "add_multiple_if_covers_more_children": True}

    if conf_parser.get_data_fetcher() == "agr_data_fetcher":
        organisms_list = conf_parser.get_agr_organisms_to_process()
    else:
        organisms_list = conf_parser.get_wb_organisms_to_process()
    for organism in organisms_list:
        logging.info("processing organism " + organism)
        sister_gene_name_id_map = {}
        sister_df = None
        species = None
        species = conf_parser.get_wb_species()
        df = WBDataFetcher(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                           release_version=conf_parser.get_release("wb_data_fetcher"),
                           species=organism, project_id=species[organism]["project_id"],
                           cache_location=conf_parser.get_cache_location(), do_relations=None,
                           go_relations=["subClassOf", "BFO:0000050"])
        if "main_sister_species" in species[organism] and species[organism]["main_sister_species"]:
            sister_df = WBDataFetcher(raw_files_source=conf_parser.get_raw_file_sources("wb_data_fetcher"),
                                      release_version=conf_parser.get_release("wb_data_fetcher"),
                                      species=species[organism]["main_sister_species"],
                                      project_id=species[species[organism]["main_sister_species"]]["project_id"],
                                      cache_location=conf_parser.get_cache_location(), do_relations=None,
                                      go_relations=["subClassOf", "BFO:0000050"])
            sister_df.load_all_data_from_file()
            for gene in sister_df.get_gene_data():
                sister_gene_name_id_map[gene.name] = gene.id
        df.load_all_data_from_file(go_terms_replacement_regex=conf_parser.get_go_rename_terms(),
                                   go_terms_exclusion_list=conf_parser.get_go_terms_exclusion_list(),
                                   do_terms_replacement_regex=None,
                                   do_terms_exclusion_list=conf_parser.get_do_terms_exclusion_list())
        desc_writer = JsonGDWriter()
        for gene in df.get_gene_data():
            logging.debug("processing gene " + gene.name)
            gene_desc = GeneDesc(gene_id=gene.id, gene_name=gene.name)
            joined_sent = []
            go_sent_generator = SentenceGenerator(
                annotations=df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                        priority_list=conf_parser.get_go_annotations_priority(),
                                                        desc_stats=gene_desc.stats),
                ontology=df.go_ontology, **go_sent_gen_common_props)
            func_sent = " and ".join([sentence.text for sentence in go_sent_generator.get_sentences(
                aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True, )])
            if func_sent:
                joined_sent.append(func_sent)
            contributes_to_func_sent = " and ".join([sentence.text for sentence in go_sent_generator.get_sentences(
                aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True,
                keep_only_best_group=True, desc_stats=gene_desc.stats, **go_sent_common_props)])
            if contributes_to_func_sent:
                joined_sent.append(contributes_to_func_sent)
            proc_sent = " and ".join([sentence.text for sentence in go_sent_generator.get_sentences(
                aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True,
                desc_stats=gene_desc.stats, **go_sent_common_props)])
            if proc_sent:
                joined_sent.append(proc_sent)
            comp_sent = " and ".join([sentence.text for sentence in go_sent_generator.get_sentences(
                aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True,
                desc_stats=gene_desc.stats, **go_sent_common_props)])
            if comp_sent:
                joined_sent.append(comp_sent)
            colocalizes_with_comp_sent = " and ".join([sentence.text for sentence in go_sent_generator.get_sentences(
                aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
                desc_stats=gene_desc.stats, keep_only_best_group=True, **go_sent_common_props)])
            if colocalizes_with_comp_sent:
                joined_sent.append(colocalizes_with_comp_sent)

            do_sentence_generator = SentenceGenerator(
                df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                            priority_list=conf_parser.get_do_annotations_priority(),
                                            desc_stats=gene_desc.stats),
                ontology=df.do_ontology, **do_sent_gen_common_prop)
            disease_sent = "; ".join([sentence.text for sentence in do_sentence_generator.get_sentences(
                aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, desc_stats=gene_desc.stats,
                **do_sent_common_props)])
            if disease_sent:
                joined_sent.append(disease_sent)

            if conf_parser.get_data_fetcher() == "wb_data_fetcher" and "main_sister_species" in species[organism] and \
                    species[organism]["main_sister_species"] and gene.name.startswith("Cbr-") and gene.name[4:] in \
                    sister_gene_name_id_map:
                sister_sentences_generator = SentenceGenerator(sister_df.get_annotations_for_gene(
                    annot_type=DataType.GO, gene_id=sister_gene_name_id_map[gene.name[4:]],
                    priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"),
                    desc_stats=gene_desc.stats), ontology=df.go_ontology, **go_sent_gen_common_props)
                sister_proc_sent = " and ".join([sentence.text for sentence in sister_sentences_generator.get_sentences(
                    aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)])
                if sister_proc_sent:
                    joined_sent.append("in " + species[species[organism]["main_sister_species"]]["name"] + ", " +
                                       gene.name[4:] + " " + sister_proc_sent)
            if len(joined_sent) > 0:
                desc = "; ".join(joined_sent) + "."
                if len(desc) > 0:
                    gene_desc.description = desc[0].upper() + desc[1:]
            else:
                gene_desc.description = None
            desc_writer.add_gene_desc(gene_desc)
        desc_writer.write(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                       organism + "_with_stats.json"), pretty=True, include_single_gene_stats=True)
        desc_writer.write(os.path.join(conf_parser.get_genedesc_output_dir(conf_parser.get_genedesc_writer()),
                                       organism + "_no_stats.json"), pretty=True, include_single_gene_stats=False)


if __name__ == '__main__':
    main()
