import argparse
import configparser
from data_fetcher import WBRawDataFetcher
from descriptions_rules import *


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="genedesc.ini", help="configuration file")
    parser.add_argument("-v", "--output-version", metavar="version_number", dest="version_number", type=str,
                        help="release version number")
    parser.add_argument("-w", "--wormbase-version", metavar="wormbase_number", dest="wormbase_number", type=str,
                        help="wormbase input files version number")

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_file)

    raw_files_source = config.get("data_fetcher", "raw_files_source")
    species = config.get("generic", "species").split(",")
    project_ids = config.get("generic", "project_ids").split(",")

    df = WBRawDataFetcher(raw_files_source=raw_files_source, release_version=args.wormbase_number, species=species[3],
                          project_id=project_ids[3])
    df.read_go_data()
    for gene in df.get_gene_data():
        go_annotations = df.get_go_annotations(gene.id, GO_ASPECT.BIOLOGICAL_PROCESS)
        sentences = generate_go_sentence([annotation.go_name for annotation in go_annotations],
                                         go_type=GO_ASPECT.BIOLOGICAL_PROCESS)
        print(sentences)


if __name__ == '__main__':
    main()