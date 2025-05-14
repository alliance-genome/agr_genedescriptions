import argparse
import logging
import time

from ontobio import Ontology

from genedescriptions.commons import DataType, Gene
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.gene_description import GeneDescription
from genedescriptions.precanned_modules import set_expression_module
from genedescriptions.descriptions_writer import DescriptionsWriter
from pipelines.alliance.alliance_data_manager import AllianceDataManager, provider_to_expression_curie_prefix
from pipelines.alliance.ateam_api_helper import get_gene_data_from_api

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.yml", help="configuration file. Default ./config.yaml")
    parser.add_argument("-L", "--log-level", dest="log_level",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default="INFO",
                        help="set the logging level")

    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
    # Set the root logger's level dynamically
    logging.getLogger(__name__).setLevel(logging.getLevelName(args.log_level))
    conf_parser = GenedescConfigParser(args.config_file)

    data_manager = AllianceDataManager(config=conf_parser)

    logger.info("Loading data providers")
    data_providers = data_manager.load_data_providers(source="db")
    data_providers.append(["HUMAN", "9606"])

    logger.info("Loading GO ontology")
    data_manager.load_ontology(ontology_type=DataType.GO, source="db")

    for data_provider, species_taxon in data_providers:
        logger.info(f"Generating gene descriptions for {data_provider}")

        json_desc_writer = DescriptionsWriter()
        if data_provider in provider_to_expression_curie_prefix:
            logger.info(f"Loading anatomy ontology data for {data_provider}")
            data_manager.load_ontology(ontology_type=DataType.EXPR, provider=data_provider, source="db")

            logger.info(f"Loading expression annotations for {data_provider}")
            data_manager.load_annotations(associations_type=DataType.EXPR, taxon_id=species_taxon,
                                          provider=data_provider, source="db")

        logger.info(f"Loading gene data for {data_provider}")
        data_manager.load_gene_data(species_taxon=species_taxon, source="db")

        logger.info(f"Generating text summaries for {data_provider}")
        for gene in data_manager.get_gene_data():
            gene_desc = GeneDescription(gene_id=gene.id,
                                        gene_name=gene.name,
                                        add_gene_name=False,
                                        config=conf_parser)
            if data_provider in provider_to_expression_curie_prefix:
                set_expression_module(df=data_manager,
                                      conf_parser=conf_parser,
                                      gene_desc=gene_desc,
                                      gene=gene)
            json_desc_writer.add_gene_desc(gene_desc)

        logger.info(f"Saving gene descriptions for {data_provider}")
        json_desc_writer.write_json(file_path=f"generated_descriptions/{data_provider}.json",
                                    include_single_gene_stats=False,
                                    data_manager=data_manager)
        json_desc_writer.write_tsv(file_path=f"generated_descriptions/{data_provider}.tsv")


if __name__ == '__main__':
    main()
