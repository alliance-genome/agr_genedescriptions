#!/usr/bin/env python3
import psycopg2 as psycopg2


def main():
    conn = psycopg2.connect("dbname='caltech_curation' user='postgres' password='' host='172.17.0.1'")
    cur = conn.cursor()
    cur.execute("select g.con_wbgene, d.con_desctext, a.con_accession, c.con_curator_hst, p.con_paper, "
                "l.con_lastupdate, per.con_person "
                "from con_wbgene g "
                "join con_desctext d ON g.joinkey = d.joinkey "
                "left outer join con_curator_hst c ON g.joinkey = c.joinkey "
                "left outer join con_paper p ON g.joinkey = p.joinkey "
                "left outer join con_accession a ON g.joinkey = a.joinkey "
                "left outer join con_lastupdate l ON g.joinkey = l.joinkey "
                "join con_desctype t ON g.joinkey = t.joinkey "
                "left outer join con_person per ON g.joinkey = per.joinkey "
                "WHERE t.con_desctype = 'Concise_description' "
                "AND g.joinkey not in ("
                "select joinkey from con_nodump) AND g.con_wbgene not in "
                "(select w.gin_wbgene from gin_dead d join gin_wbgene w ON d.joinkey = w.joinkey)")
    rows = cur.fetchall()
    genedesc = {}
    for row in rows:
        if row[0] in genedesc and row[3]:
            genedesc[row[0]][2].add(row[3])
        else:
            genedesc[row[0]] = [row[1].replace("\n", ""), row[2], set([row[3]] if row[3] else []), row[4], row[5], row[6]]

    for gene_id, gene_props in genedesc.items():
        desc_text = gene_props[0]
        print("Gene : \"" + gene_id + "\"")
        if not gene_props[1] and not gene_props[2] and not len(gene_props[4]) > 0 and not gene_props[5]:
            print("Concise_description", "\"" + desc_text + "\"", sep="\t")
        if gene_props[1]:
            for accession in gene_props[1].split(", "):
                accession_arr = accession.split(":")
                print("Concise_description", "\"" + desc_text + "\"", "Accession_evidence", "\"" + accession_arr[0] +
                      "\" \"" + accession_arr[1] + "\"", sep="\t")
        if gene_props[2]:
            for person in gene_props[2]:
                print("Concise_description", "\"" + desc_text + "\"", "Curator_confirmed", "\"" + person + "\"",
                      sep="\t")
        if gene_props[3]:
            for paper in gene_props[3].split(","):
                print("Concise_description", "\"" + desc_text + "\"", "Paper_evidence", paper, sep="\t")
        if gene_props[4]:
            print("Concise_description", "\"" + desc_text + "\"", "Date_last_updated", "\"" + gene_props[4].split(" ")[0] +
                  "\"", sep="\t")
        if gene_props[5]:
            for person in gene_props[5].split(","):
                print("Concise_description", "\"" + desc_text + "\"", "Person_evidence", person, sep="\t")
        print()


if __name__ == '__main__':
    main()
