#!/usr/bin/env python3
import psycopg2 as psycopg2


def main():
    conn = psycopg2.connect("dbname='testdb' user='acedb' password='' host='tazendra.caltech.edu'")
    cur = conn.cursor()
    cur.execute("select g.con_wbgene, d.con_desctext, a.con_accession, c.con_curator, p.con_paper, l.con_lastupdate, "
                "per.con_person "
                "from con_wbgene g "
                "join con_desctext d ON g.joinkey = d.joinkey "
                "left outer join con_curator c ON g.joinkey = c.joinkey "
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
    for row in rows:
        desc_text = row[1].replace("\n", "")
        print("Gene : \"" + row[0] + "\"")
        if not row[2] and not row[3] and not row[4] and not row[6]:
            print("Concise_description", "\"" + desc_text + "\"", sep="\t")
        if row[2]:
            for accession in row[2].split(", "):
                accession_arr = accession.split(":")
                print("Concise_description", "\"" + desc_text + "\"", "Accession_evidence", "\"" + accession_arr[0] +
                      "\" \"" + accession_arr[1] + "\"", sep="\t")
        if row[3]:
            for person in row[3].split(", "):
                print("Concise_description", "\"" + desc_text + "\"", "Curator_confirmed", "\"" + person + "\"",
                      sep="\t")
        if row[4]:
            for paper in row[4].split(","):
                print("Concise_description", "\"" + desc_text + "\"", "Paper_evidence", paper, sep="\t")
        if row[5]:
            print("Concise_description", "\"" + desc_text + "\"", "Date_last_updated", "\"" + row[5].split(" ")[0] +
                  "\"", sep="\t")
        if row[6]:
            for person in row[6].split(","):
                print("Concise_description", "\"" + desc_text + "\"", "Person_evidence", person, sep="\t")
        print()


if __name__ == '__main__':
    main()
