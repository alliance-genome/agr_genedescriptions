import argparse
import configparser
import urllib
import json

from data_fetcher import DataFetcher


def generate_go_sentences(go_type: str, gene, sentence_header):
    """"""
    if go_type in gene:
        functions = []
        for function in gene[go_type]:
            functions.append(function)
        if len(functions) > 1:
            return sentence_header + " " + ", ".join(functions[0:-1]) + ", and " + functions[len(functions) - 1]
        else:
            return sentence_header + " " + functions[0]


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

    df = DataFetcher(raw_files_source=raw_files_source, release_version=args.wormbase_number, species=species[0],
                     project_id=project_ids[0])
    for geneid in df.get_gene_data():
        print(geneid)

    with urllib.request.urlopen(
            "http://www.alliancegenome.org/api/search?category=gene&limit=100&offset=0&species=" + species) as url:
        data = json.loads(url.read().decode())
        for gene in data["results"]:
            # header
            print(gene["id"], gene["symbol"])
            # description
            sentences = []
            mol_func = generate_go_sentences("gene_molecular_function", gene, "encodes product that exhibits")
            if mol_func:
                sentences.append(mol_func)
            bio_proc = generate_go_sentences("gene_biological_process", gene, "is involved in")
            if bio_proc:
                sentences.append(bio_proc)
            cell_comp = generate_go_sentences("gene_cellular_component", gene, "localizes to")
            if cell_comp:
                sentences.append(cell_comp)
            final_sentence = "; ".join(sentences).capitalize()
            if len(final_sentence) == 0:
                final_sentence = "No Description"
            final_sentence += "."
            print(final_sentence)
            print()


if __name__ == '__main__':
    main()