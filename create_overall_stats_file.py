#!/usr/bin/env python3
import argparse
import os
import json


def main():
    parser = argparse.ArgumentParser(description="Combine automated description statistics from json files and manual "
                                                 "descriptions stats from concise descriptions ace file")
    parser.add_argument("-i", "--input-dir", metavar="input_dir", dest="input_dir", type=str,
                        default="./", help="working directory where input json files are located")
    parser.add_argument("-m", "--manual-desc", metavar="manual_desc", dest="manual_desc", type=str,
                        help="path to the manual descriptions ace file")
    args = parser.parse_args()

    genes_with_automated_desc = set()
    genes_with_manual_desc = set([line.split(" : ")[1].replace("\"", "") for line in open(args.manual_desc) if
                                  line.startswith("Gene : ")])
    num_orthology_sentences = 0
    num_go_process_sentences = 0
    num_go_function_sentences = 0
    num_go_component_sentences = 0
    num_tissue_expression_sentences = 0
    num_gene_expression_cluster_sentences = 0
    num_molecule_expression_cluster_sentences = 0
    num_anatomy_expression_cluster_sentences = 0
    num_disease_experimental_sentences = 0
    num_disease_biomarker_sentences = 0
    num_disease_orthology_sentences = 0
    num_protein_domain_sentences = 0
    num_human_gene_function_sentences = 0
    num_sister_species_sentences = 0

    for json_file_path in os.listdir(args.input_dir):
        if json_file_path.endswith(".json"):
            with open(json_file_path) as json_file:
                json_data = json.loads(json_file)
                genes_with_automated_desc.update(set([desc["gene_id"] for desc in json_data["data"] if
                                                      desc["description"] is not None]))
                num_orthology_sentences += len([desc for desc in json_data["data"] if desc["orthology_description"]
                                                is not None])
                num_go_process_sentences += len([desc for desc in json_data["data"] if desc["go_process_description"]
                                                is not None])
                num_go_function_sentences += len([desc for desc in json_data["data"] if desc["go_function_description"]
                                                  is not None])
                num_go_component_sentences += len([desc for desc in json_data["data"] if
                                                   desc["go_component_description"] is not None])
                num_tissue_expression_sentences += len([desc for desc in json_data["data"] if
                                                        desc["tissue_expression_description"] is not None])
                num_gene_expression_cluster_sentences += len([desc for desc in json_data["data"] if
                                                              desc["gene_expression_cluster_description"] is not None])
                num_molecule_expression_cluster_sentences += len([desc for desc in json_data["data"] if
                                                                  desc["molecule_expression_cluster_description"] is not
                                                                  None])
                num_anatomy_expression_cluster_sentences += len([desc for desc in json_data["data"] if
                                                                 desc["anatomy_expression_cluster_description"] is not
                                                                 None])
                num_disease_experimental_sentences += len([desc for desc in json_data["data"] if
                                                           desc["do_experimental_description"] is not None])
                num_disease_biomarker_sentences += len([desc for desc in json_data["data"] if
                                                        desc["do_biomarker_description"] is not None])
                num_disease_orthology_sentences += len([desc for desc in json_data["data"] if
                                                        desc["do_orthology_description"] is not None])
                num_protein_domain_sentences += len([desc for desc in json_data["data"] if
                                                     desc["protein_domain_description"] is not None])
                num_human_gene_function_sentences += len([desc for desc in json_data["data"] if
                                                          desc["human_gene_function_description"] is not None])
                num_sister_species_sentences += len([desc for desc in json_data["data"] if
                                                     desc["sister_species_description"] is not None])

    genes_with_both_desc = len(genes_with_automated_desc.intersection(genes_with_manual_desc))

    # TODO add total number of automated and concise descriptions
    # TODO add number of genes with both concise and automated




if __name__ == '__main__':
    main()
