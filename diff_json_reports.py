#!/usr/bin/env python3

import argparse
import json


def main():
    parser = argparse.ArgumentParser(description="Create diff file between reports in json format")
    parser.add_argument('reports', metavar='report_file', type=str, nargs=2, help='report file to diff')

    args = parser.parse_args()

    old_report = args.reports[0]
    new_report = args.reports[1]

    old_gene_descs = {}
    diff_descs = []

    with open(old_report, 'r') as f:
        old_report_data = json.load(f)["data"]
    for genedesc in old_report_data:
        old_gene_descs[genedesc["gene_id"]] = genedesc

    with open(new_report, 'r') as f:
        new_report_data = json.load(f)["data"]

    for genedesc in new_report_data:
        if genedesc["gene_id"] not in old_gene_descs:
            new_genedesc = {"gene_id": genedesc["gene_id"],
                            "gene_name": genedesc["gene_name"],
                            "old_description": None,
                            "new_description": genedesc["description"],
                            "new_number_go_annotations":
                                genedesc["stats"]["total_number_go_annotations"],
                            "old_number_go_annotations": 0,
                            "new_set_initial_go_ids_f":
                                genedesc["stats"]["set_initial_go_ids_f"],
                            "old_set_initial_go_ids_f": 0,
                            "new_set_initial_go_ids_p":
                                genedesc["stats"]["set_initial_go_ids_p"],
                            "old_set_initial_go_ids_p": 0,
                            "new_set_initial_go_ids_c":
                                genedesc["stats"]["set_initial_go_ids_c"],
                            "old_set_initial_go_ids_c": 0,
                            "new_set_final_go_ids_f":
                                genedesc["stats"]["set_final_go_ids_f"],
                            "old_set_final_go_ids_f": 0,
                            "new_set_final_go_ids_p":
                                genedesc["stats"]["set_final_go_ids_p"],
                            "old_set_final_go_ids_p": 0,
                            "new_set_final_go_ids_c":
                                genedesc["stats"]["set_final_go_ids_c"],
                            "old_set_final_go_ids_c": 0,
                            "new_number_do_annotations":
                                genedesc["stats"]["total_number_do_annotations"],
                            "old_number_do_annotations": 0,
                            "new_number_initial_do_terms":
                                genedesc["stats"]["number_initial_do_terms"],
                            "old_number_initial_do_terms": 0,
                            "new_number_final_do_terms":
                                genedesc["stats"]["number_final_do_terms"],
                            "old_number_final_do_terms": 0,
                            "new_number_final_do_term_covering_multiple_initial_do_terms":
                                genedesc["stats"]["number_final_do_term_covering_multiple_initial_do_terms"],
                            "old_number_final_do_term_covering_multiple_initial_do_terms": 0}
            diff_descs.append(new_genedesc)
        elif genedesc["description"] != old_gene_descs[genedesc["gene_id"]]["description"]:
            new_genedesc = {"gene_id": genedesc["gene_id"],
                            "gene_name": genedesc["gene_name"],
                            "old_description":  old_gene_descs[genedesc["gene_id"]]["description"],
                            "new_description": genedesc["description"],
                            "new_number_go_annotations": genedesc["stats"]["total_number_go_annotations"],
                            "old_number_go_annotations": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "total_number_go_annotations"],
                            "new_set_initial_go_ids_f": genedesc["stats"]["set_initial_go_ids_f"],
                            "old_set_initial_go_ids_f": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_initial_go_ids_f"],
                            "new_set_initial_go_ids_p": genedesc["stats"]["set_initial_go_ids_p"],
                            "old_set_initial_go_ids_p": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_initial_go_ids_p"],
                            "new_set_initial_go_ids_c": genedesc["stats"]["set_initial_go_ids_c"],
                            "old_set_initial_go_ids_c": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_initial_go_ids_c"],
                            "new_set_final_go_ids_f": genedesc["stats"]["set_final_go_ids_f"],
                            "old_set_final_go_ids_f": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_final_go_ids_f"],
                            "new_set_final_go_ids_p": genedesc["stats"]["set_final_go_ids_p"],
                            "old_set_final_go_ids_p": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_final_go_ids_p"],
                            "new_set_final_go_ids_c": genedesc["stats"]["set_final_go_ids_c"],
                            "old_set_final_go_ids_c": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "set_final_go_ids_c"],
                            "new_number_do_annotations": genedesc["stats"]["total_number_do_annotations"],
                            "old_number_do_annotations": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "total_number_do_annotations"],
                            "new_number_initial_do_terms": genedesc["stats"]["number_initial_do_terms"],
                            "old_number_initial_do_terms": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "number_initial_do_terms"],
                            "new_number_final_do_terms": genedesc["stats"]["number_final_do_terms"],
                            "old_number_final_do_terms": old_gene_descs[genedesc["gene_id"]]["stats"][
                                "number_final_do_terms"],
                            "new_number_final_do_term_covering_multiple_initial_do_terms":
                                genedesc["stats"]["number_final_do_term_covering_multiple_initial_do_terms"],
                            "old_number_final_do_term_covering_multiple_initial_do_terms":
                                old_gene_descs[genedesc["gene_id"]]["stats"][
                                    "number_final_do_term_covering_multiple_initial_do_terms"]}
            diff_descs.append(new_genedesc)

    print(json.dumps(diff_descs, indent=True))


if __name__ == '__main__':
    main()
