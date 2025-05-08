import argparse
import logging
import time

from genedescriptions.commons import DataType, Gene
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.gene_description import GeneDescription
from genedescriptions.precanned_modules import set_expression_module
from genedescriptions.descriptions_writer import DescriptionsWriter
from pipelines.alliance.alliance_data_manager import AllianceDataManager


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
    logger.info("Loading anatomy ontology data")
    logger.info("Measuring time to load ontology from API")
    start_time_api = time.time()
    data_manager.load_ontology_from_ateam_api(ontology_type=DataType.EXPR, provider="WB")
    end_time_api = time.time()
    time_api = end_time_api - start_time_api
    logger.info(f"Time to load ontology from API: {time_api:.2f} seconds")

    logger.info("Measuring time to load ontology from persistent store")
    start_time_store = time.time()
    data_manager.load_ontology_from_persistent_store(ontology_type=DataType.EXPR, provider="WB")
    end_time_store = time.time()
    time_store = end_time_store - start_time_store
    logger.info(f"Time to load ontology from persistent store: {time_store:.2f} seconds")

    logger.info("Loading gene data")
    data_manager.load_gene_data_from_persistent_store(provider="WB")
    logger.info("Loading expression annotations")
    data_manager.load_annotations_from_persistent_store(associations_type=DataType.EXPR,
                                                        taxon_id="NCBITaxon:6239", provider="WB")
    json_desc_writer = DescriptionsWriter()
    for gene in data_manager.get_gene_data():
        gene_desc = GeneDescription(gene_id=gene.id,
                                    gene_name=gene.name,
                                    add_gene_name=False,
                                    config=conf_parser)
        set_expression_module(df=data_manager,
                              conf_parser=conf_parser,
                              gene_desc=gene_desc,
                              gene=gene)
        json_desc_writer.add_gene_desc(gene_desc)

    json_desc_writer.write_json(file_path="wormbase.json", include_single_gene_stats=False, data_manager=data_manager)


if __name__ == '__main__':
    main()
