#!/usr/bin/env bash

python3 wormbase_pipeline.py -c config_wb.yml -L DEBUG -o json txt tsv ace -t $1
cat generated_descriptions/*.ace > automatedesc.ace
rm generated_descriptions/*.ace
mv automatedesc.ace generated_descriptions/
python3 manualdesc_postgres2ace.py > generated_descriptions/concisedesc.ace
python3 create_overall_stats_file.py -i generated_descriptions -m generated_descriptions/concisedesc.ace > generated_descriptions/__number_of_concise_descriptions.txt
