import urllib
import json


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
    species = "Mus+musculus"
    limit = 100
    offset = 0
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