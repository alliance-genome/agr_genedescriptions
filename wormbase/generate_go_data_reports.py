import argparse
import sys

from ontobio import OntologyFactory, AssociationSetFactory
from ontobio.io.assocparser import AssocParserConfig
from ontobio.io.gafparser import GafParser


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Load GO ontology and GAF files to get and print error reports.")
    parser.add_argument("-o", "--ontology_file", help="Path to the GO ontology file to be loaded.",
                        required=True)
    parser.add_argument("-g", "--gaf_file", help="Path to the GAF file to be loaded.", required=True)
    args = parser.parse_args()

    try:
        # Load GO ontology file
        print(f"Loading ontology file: {args.ontology_file}")
        go_ontology = OntologyFactory().create(args.ontology_file)

        # Load GAF file
        print(f"Loading GAF file: {args.gaf_file}")
        assoc_config = AssocParserConfig(remove_double_prefixes=True, paint=True)
        gaf_parser = GafParser(config=assoc_config)
        assocs = AssociationSetFactory().create_from_assocs(
            assocs=gaf_parser.parse(file=args.gaf_file, skipheader=True),
            ontology=go_ontology)

        # Retrieve error reports
        print("Retrieving and printing error reports...")
        errors = gaf_parser.report.messages
        for error in errors:
            print(f"ERROR: {error}")

        print("Processing complete.")

    except Exception as e:
        print(f"An error occurred during processing: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
