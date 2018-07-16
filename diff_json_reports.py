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
        if genedesc["description"] != old_gene_descs[genedesc["gene_id"]]["description"]:
            diff_descs.append(genedesc)

    print(json.dumps(diff_descs, indent=True))


if __name__ == '__main__':
    main()
