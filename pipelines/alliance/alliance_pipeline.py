import argparse
import logging

from genedescriptions.commons import DataType
from genedescriptions.config_parser import GenedescConfigParser
from pipelines.alliance.alliance_data_manager import AllianceDataManager


logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Generate gene descriptions for wormbase")
    parser.add_argument("-c", "--config-file", metavar="config_file", dest="config_file", type=str,
                        default="config.yml", help="configuration file. Default ./config.yaml")
    parser.add_argument("-L", "--log-level", dest="log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                                                        'CRITICAL'], help="set the logging level")

    args = parser.parse_args()
    conf_parser = GenedescConfigParser(args.config_file)

    data_manager = AllianceDataManager(config=conf_parser)
    data_manager.load_gene_data_from_persistent_store(provider="WB")
    data_manager.load_ontology_from_persistent_store(ontology_type=DataType.EXPR, provider="WB")
    data_manager.load_annotations_from_persistent_store(associations_type=DataType.EXPR,
                                                        taxon_id="NCBITaxon:6239", provider="WB")




if __name__ == '__main__':
    main()
