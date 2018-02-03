import argparse
import configparser
from data_fetcher import WBRawDataFetcher
from descriptions_rules import *


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="genedesc.ini", help="configuration file")
    parser.add_argument("-C", "--use-cache", dest="use_cache", action="store_true",
                        help="Use cached source files from cache_location specified in config file. Download them from "
                             "raw_file_source (configured in config file) if not yet cached")
    parser.add_argument("-v", "--output-version", metavar="version_number", dest="version_number", type=str,
                        help="release version number")
    parser.add_argument("-w", "--wormbase-version", metavar="wormbase_number", dest="wormbase_number", type=str,
                        help="wormbase input files version number")

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config_file)

    raw_files_source = config.get("data_fetcher", "raw_files_source")
    cache_location = config.get("data_fetcher", "cache_location")
    species = config.get("generic", "species").split(",")
    project_ids = config.get("generic", "project_ids").split(",")

    df = WBRawDataFetcher(raw_files_source=raw_files_source, release_version=args.wormbase_number, species=species[3],
                          project_id=project_ids[3], cache_location=cache_location, use_cache=True)
    df.load_go_data()
    for gene in df.get_gene_data():
        print(gene.id, gene.name)
        sentences = generate_go_sentence(df.get_go_annotations(gene.id))
        if sentences:
            joined_sent = []
            proc_sent = "; ".join([sent.text for sent in sentences.get_sentences(GO_ASPECT.BIOLOGICAL_PROCESS)])
            if proc_sent:
                joined_sent.append(proc_sent.capitalize())
            func_sent = "; ".join([sent.text for sent in sentences.get_sentences(GO_ASPECT.MOLECULAR_FUNCTION)])
            if func_sent:
                joined_sent.append(func_sent.capitalize())
            comp_sent = "; ".join([sent.text for sent in sentences.get_sentences(GO_ASPECT.CELLULAR_COMPONENT)])
            if comp_sent:
                joined_sent.append(comp_sent.capitalize())

            print(". ".join(joined_sent) + ".")
        else:
            print("No description available")
        print()


if __name__ == '__main__':
    main()
