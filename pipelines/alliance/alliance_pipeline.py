import argparse
import concurrent.futures
import logging
import time

from genedescriptions.commons import DataType
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.descriptions_writer import DescriptionsWriter
from genedescriptions.gene_description import GeneDescription
from genedescriptions.precanned_modules import set_expression_module, set_gene_ontology_module
from pipelines.alliance.alliance_data_manager import AllianceDataManager, provider_to_expression_curie_prefix

logger = logging.getLogger(__name__)

DATA_SOURCE = "db"


def load_all_data_for_provider(data_manager: AllianceDataManager, data_provider: str, species_taxon: str):
    logger.info(f"Loading GAF file for {data_provider}")
    data_manager.load_annotations(associations_type=DataType.GO, taxon_id=species_taxon, provider=data_provider,
                                  source=DATA_SOURCE)
    if data_provider in provider_to_expression_curie_prefix:
        logger.info(f"Loading anatomy ontology data for {data_provider}")
        data_manager.load_ontology(ontology_type=DataType.EXPR, provider=data_provider, source="db")

        logger.info(f"Loading expression annotations for {data_provider}")
        data_manager.load_annotations(associations_type=DataType.EXPR, taxon_id=species_taxon,
                                      provider=data_provider, source=DATA_SOURCE)

    logger.info(f"Loading gene data for {data_provider}")
    data_manager.load_gene_data(species_taxon=species_taxon, source=DATA_SOURCE)


def generate_gene_descriptions(data_manager: AllianceDataManager, data_provider: str,
                               conf_parser: GenedescConfigParser, json_desc_writer: DescriptionsWriter):
    for gene in data_manager.get_gene_data():
        gene_desc = GeneDescription(gene_id=gene.id,
                                    gene_name=gene.name,
                                    add_gene_name=False,
                                    config=conf_parser)
        set_gene_ontology_module(dm=data_manager, conf_parser=conf_parser, gene_desc=gene_desc, gene=gene)
        if data_provider in provider_to_expression_curie_prefix:
            set_expression_module(df=data_manager,
                                  conf_parser=conf_parser,
                                  gene_desc=gene_desc,
                                  gene=gene)
        json_desc_writer.add_gene_desc(gene_desc)


def save_gene_descriptions(data_manager: AllianceDataManager, json_desc_writer: DescriptionsWriter, data_provider: str):
    json_desc_writer.write_json(file_path=f"generated_descriptions/{data_provider}.json",
                                include_single_gene_stats=False,
                                data_manager=data_manager)
    json_desc_writer.write_tsv(file_path=f"generated_descriptions/{data_provider}.tsv")


def process_provider(data_provider, species_taxon, data_manager, conf_parser):
    logger.info(f"Processing provider: {data_provider}")
    json_desc_writer = DescriptionsWriter()

    logger.info(f"Loading all data for {data_provider}")
    load_all_data_for_provider(data_manager, data_provider, species_taxon)

    logger.info(f"Generating text summaries for {data_provider}")
    generate_gene_descriptions(data_manager, data_provider, conf_parser, json_desc_writer)

    logger.info(f"Saving gene descriptions for {data_provider}")
    save_gene_descriptions(data_manager, json_desc_writer, data_provider)


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.yml", help="configuration file. Default ./config.yaml")
    parser.add_argument("-L", "--log-level", dest="log_level",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default="INFO",
                        help="set the logging level")
    parser.add_argument("--parallel", dest="parallel", action="store_true",
                        help="Run providers in parallel")
    parser.add_argument("--max-workers", dest="max_workers", type=int, default=None,
                        help="Maximum number of parallel executions (default: number of CPUs)")

    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
    logging.getLogger(__name__).setLevel(logging.getLevelName(args.log_level))

    start_time = time.time()
    conf_parser = GenedescConfigParser(args.config_file)

    data_manager = AllianceDataManager(config=conf_parser)

    logger.info("Loading data providers")
    data_providers = data_manager.load_data_providers(source=DATA_SOURCE)
    data_providers.append(["HUMAN", "9606"])

    logger.info("Loading GO ontology")
    data_manager.load_ontology(ontology_type=DataType.GO, source=DATA_SOURCE)

    if args.parallel:
        logger.info("Processing data providers in parallel")
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers) as executor:
            futures = [
                executor.submit(process_provider, data_provider, species_taxon, data_manager, conf_parser)
                for data_provider, species_taxon in data_providers
            ]
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Error processing data provider: {e}")
    else:
        logger.info("Processing data providers sequentially")
        for data_provider, species_taxon in data_providers:
            try:
                process_provider(data_provider, species_taxon, data_manager, conf_parser)
            except Exception as e:
                logger.error(f"Error processing data provider: {e}")

    elapsed_time = time.time() - start_time
    formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    logger.info(f"All data providers processed successfully in {formatted_time}")


if __name__ == '__main__':
    main()
