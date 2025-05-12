#!/usr/bin/env python3
import argparse


def main():
    parser = argparse.ArgumentParser(description="Return all gene entries in ace file that are present in the provided "
                                                 "WB gene list file")
    parser.add_argument("-a", "--ace-file", metavar="ace_file", dest="ace_file", type=str,
                        help="a valid ace file")
    parser.add_argument("-g", "--genes-file", dest="genes_file", metavar="genes_file", type=str,
                        help="file containing gene list")
    args = parser.parse_args()

    valid_gene_ids = set([line.strip().split(",")[1] for line in open(args.genes_file) if
                          line.strip().split(",")[4] == "Live"])

    print_lines = []
    gene_id = None
    for line in open(args.ace_file):
        if line == "\n":
            if gene_id in valid_gene_ids and print_lines:
                print()
                for ln in print_lines:
                    print(ln)
            print_lines = []
            gene_id = None
        else:
            linearr = line.strip().split(" ")
            if len(linearr) == 3 and linearr[0] == "Gene":
                gene_id = linearr[2][1:-1]
            print_lines.append(line.strip())


if __name__ == '__main__':
    main()
