#!/usr/bin/env bash

export PYTHONPATH=$PYTHONPATH:/usr/src/gene_descriptions

python3 /usr/src/gene_descrptions/wormbase_pipeline.py -c /usr/src/gene_descrptions/wormbase/config_wb.yml -L DEBUG -o json txt tsv ace -t $1
cat /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/*.ace > automatedesc.ace
rm /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/*.ace
today_folder=$(date '+%Y%m%d')
mkdir /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder}
mv automatedesc.ace /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder}
python3 /usr/src/gene_descrptions/manualdesc_postgres2ace.py > /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}_manualdesc.ace
python3 /usr/src/gene_descrptions/create_overall_stats_file.py -i /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder} -m /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}_manualdesc.ace > /usr/caltech_curation_files/pub/gene_descrpitons/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}__number_of_concise_descriptions.txt
