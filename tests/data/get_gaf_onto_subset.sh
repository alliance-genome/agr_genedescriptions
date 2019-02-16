#!/bin/bash

if [ "$#" -ne 1 ]
then
    echo "missing argument: folder containing original gaf files"
    exit 1
fi

genes='HTT\|Htt\|htt\|CFTR\|Cftr\|cftr\|YTA12\|tefu\|ATM\|ATM1\|Atm\|MCM1\|Srf\|SRF\|srfb\|REB1\|daf-16\|aap-1\|abl-1\|pat-4\|Pten\|PTEN\|Ahr\|AHR\|CDC28\|CDC48\|CDC5\|SGS1\|LAS17\|RAD27\|MLH1\|RAD53\|CCC2\|COQ6\|TRS20\|PXA1\|DPB3\|lsd-1\|let-60\|MPZ\|ntn1b\|let-60\|alg-1\|UBP6\|GPM2'

for assoc_file in $1/*
do
    grep ${genes} "${assoc_file}" > $(basename "${assoc_file}").partial
done

read -r -d '' HEADER <<- EOM
[URI of the OWL(RDF/XML) output file]
http://purl.obolibrary.org/obo/go_gd_tests.owl

[Source ontology]
#comment here
GO

[Low level source term URIs]
EOM

read -r -d '' FOOTER <<- EOM
[Top level source term URIs and target direct superclass URIs]
http://purl.obolibrary.org/obo/GO_0008150 #biological_process
http://purl.obolibrary.org/obo/GO_0003674 #molecular_function
http://purl.obolibrary.org/obo/GO_0005575 #cellular_component

[Source term retrieval setting]
includeAllIntermediates

[Source annotation URIs]
includeAllAxioms
EOM

echo "${HEADER}" > ontofox_script.txt

grep -v '^!.*$' $1/* | grep ${genes} | cut -f5 | sort | uniq | sed 's/:/_/' | awk '{print "http://purl.obolibrary.org/obo/"$1}' >> ontofox_script.txt

echo "http://purl.obolibrary.org/obo/GO_0043055" >> ontofox_script.txt
echo "http://purl.obolibrary.org/obo/GO_0061065" >> ontofox_script.txt
echo "http://purl.obolibrary.org/obo/GO_0043054" >> ontofox_script.txt
echo "http://purl.obolibrary.org/obo/GO_0043053" >> ontofox_script.txt
echo "http://purl.obolibrary.org/obo/GO_0006096" >> ontofox_script.txt

echo >> ontofox_script.txt

echo "${FOOTER}" >> ontofox_script.txt

curl -s -F file=@"./ontofox_script.txt" -o go_gd_test.owl http://ontofox.hegroup.org/service.php
rm ontofox_script.txt
robot convert --input go_gd_test.owl --format obo --output go_gd_test.obo
rm go_gd_test.owl
