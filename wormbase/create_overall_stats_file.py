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
                        help="path to the manual descriptions ace file", required=True)
    args = parser.parse_args()

    genes_with_automated_desc = set()
    genes_with_manual_desc = set([line.split(" : ")[1].replace("\"", "").strip() for line in open(args.manual_desc) if
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

    for json_file_path in sorted(os.listdir(args.input_dir)):
        if json_file_path.endswith(".json"):
            with open(os.path.join(args.input_dir, json_file_path)) as json_file:
                json_data = json.load(json_file)
                genes_with_non_null_descriptions = set([desc["gene_id"].replace("WB:", "") for desc in json_data["data"]
                                                        if desc["description"] is not None])
                partial_num_orthology_sentences = len([desc for desc in json_data["data"] if
                                                       desc["orthology_description"] is not None])
                partial_num_go_process_sentences = len([desc for desc in json_data["data"] if
                                                        desc["go_process_description"] is not None])
                partial_num_go_function_sentences = len([desc for desc in json_data["data"] if
                                                         desc["go_function_description"] is not None])
                partial_num_go_component_sentences = len([desc for desc in json_data["data"] if
                                                          desc["go_component_description"] is not None])
                partial_num_tissue_expression_sentences = len([desc for desc in json_data["data"] if
                                                               desc["tissue_expression_description"] is not None])
                partial_num_gene_expression_cluster_sentences = len([desc for desc in json_data["data"] if
                                                                     desc["gene_expression_cluster_description"] is not
                                                                     None])
                partial_num_molecule_expression_cluster_sentences = len([desc for desc in json_data["data"] if
                                                                         desc["molecule_expression_cluster_description"]
                                                                         is not None])
                partial_num_anatomy_expression_cluster_sentences = len([desc for desc in json_data["data"] if
                                                                        desc["anatomy_expression_cluster_description"]
                                                                        is not None])
                partial_num_disease_experimental_sentences = len([desc for desc in json_data["data"] if
                                                                  desc["do_experimental_description"] is not None])
                partial_num_disease_biomarker_sentences = len([desc for desc in json_data["data"] if
                                                               desc["do_biomarker_description"] is not None])
                partial_num_disease_orthology_sentences = len([desc for desc in json_data["data"] if
                                                               desc["do_orthology_description"] is not None])
                partial_num_protein_domain_sentences = len([desc for desc in json_data["data"] if
                                                            desc["protein_domain_description"] is not None])
                partial_num_human_gene_function_sentences = len([desc for desc in json_data["data"] if
                                                                 desc["human_gene_function_description"] is not None])
                partial_num_sister_species_sentences = len([desc for desc in json_data["data"] if
                                                            desc["sister_species_description"] is not None])
                genes_with_automated_desc.update(genes_with_non_null_descriptions)
                num_orthology_sentences += partial_num_orthology_sentences
                num_go_process_sentences += partial_num_go_process_sentences
                num_go_function_sentences += partial_num_go_function_sentences
                num_go_component_sentences += partial_num_go_component_sentences
                num_tissue_expression_sentences += partial_num_tissue_expression_sentences
                num_gene_expression_cluster_sentences += partial_num_gene_expression_cluster_sentences
                num_molecule_expression_cluster_sentences += partial_num_molecule_expression_cluster_sentences
                num_anatomy_expression_cluster_sentences += partial_num_anatomy_expression_cluster_sentences
                num_disease_experimental_sentences += partial_num_disease_experimental_sentences
                num_disease_biomarker_sentences += partial_num_disease_biomarker_sentences
                num_disease_orthology_sentences += partial_num_disease_orthology_sentences
                num_protein_domain_sentences += partial_num_protein_domain_sentences
                num_human_gene_function_sentences += partial_num_human_gene_function_sentences
                num_sister_species_sentences += partial_num_sister_species_sentences

                print(json_data["overall_properties"]["species"])
                print(str(len(genes_with_non_null_descriptions)) + " individual gene descriptions")
                print(str(len(genes_with_non_null_descriptions.intersection(genes_with_manual_desc))) +
                      " genes have manual descriptions")
                print(str(len(genes_with_non_null_descriptions.difference(genes_with_manual_desc))) +
                      " genes have only automated descriptions")
                print(str(partial_num_orthology_sentences) + " orthology sentences")
                print(str(partial_num_go_process_sentences) + " gene ontology process sentences")
                print(str(partial_num_go_function_sentences) + " gene ontology molecular function sentences")
                print(str(partial_num_go_component_sentences) + " gene ontology cellular component sentences")
                print(str(partial_num_tissue_expression_sentences) + " tissue expression sentences")
                print(str(partial_num_gene_expression_cluster_sentences) + " gene regulation expression cluster sentences")
                print(str(partial_num_molecule_expression_cluster_sentences) + " molecule regulation expression cluster sentences")
                print(str(partial_num_anatomy_expression_cluster_sentences) + " gene expression (anatomy) cluster sentences")
                print(str(partial_num_disease_experimental_sentences) + " disease sentences based on experimental data")
                print(str(partial_num_disease_biomarker_sentences) + " disease sentences based on biomarker data")
                print(str(partial_num_disease_orthology_sentences) + " disease sentences based on orthology data")
                print(str(partial_num_protein_domain_sentences) + " protein domain sentences")
                print(str(partial_num_human_gene_function_sentences) + " human gene GO molecular function sentences")
                print(str(partial_num_sister_species_sentences) + " elegans process sentences in non-elegans species")
                print()

    print()
    print("Total number of genes with either automated or manual description: " + str(len(genes_with_automated_desc.union(
        genes_with_manual_desc))))
    print("Total number of genes with only automated description: " + str(len(
        genes_with_automated_desc.difference(genes_with_manual_desc))))
    print("Total number of genes with automated description (with or without manual description): " + str(len(
        genes_with_automated_desc)))
    print("Total number of genes with only manual description: " + str(len(
        genes_with_manual_desc.difference(genes_with_automated_desc))))
    print("Total number of genes with manual description (with or without automated description): " + str(len(
        genes_with_manual_desc)))
    print("Total number of genes with both manual and automated description: " + str(len(
        genes_with_manual_desc.intersection(genes_with_automated_desc))))


if __name__ == '__main__':
    main()
