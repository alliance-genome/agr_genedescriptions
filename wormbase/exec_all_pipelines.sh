export PYTHONPATH=$PYTHONPATH:/usr/src/gene_descriptions

today_folder=$(date '+%Y%m%d')
mkdir -p /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}
python3 /usr/src/gene_descriptions/wormbase/wormbase_pipeline.py -c /usr/src/gene_descriptions/wormbase/config_wb.yml -L DEBUG -o json txt tsv ace -t ${TEXTPRESSO_TOKEN}
cat /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/*.ace > /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}_automatedesc.ace
rm /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/*.ace
mv /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/*.{json,txt,tsv} /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}
python3 /usr/src/gene_descriptions/wormbase/manualdesc_postgres2ace.py > /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}_manualdesc.ace
python3 /usr/src/gene_descriptions/wormbase/create_overall_stats_file.py -i /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder} -m /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}_manualdesc.ace > /usr/caltech_curation_files/pub/gene_descriptions/${WB_RELEASE}/pre-release/${today_folder}/${today_folder}__number_of_concise_descriptions.txt
