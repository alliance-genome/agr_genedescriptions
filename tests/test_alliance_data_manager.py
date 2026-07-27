import logging
import os
import sys
import unittest
from unittest.mock import MagicMock

from ontobio import Ontology

from genedescriptions.config_parser import GenedescConfigParser

sys.path.insert(0, os.path.join(os.path.split(__file__)[0], os.path.pardir, "pipelines", "alliance"))

from alliance_data_manager import AllianceDataManager, GENE_SO_TERMS  # noqa: E402

logger = logging.getLogger("Alliance DataManager tests")


class TestAllianceDataManagerObsolete(unittest.TestCase):
    """Tests for obsolete ontology term handling in the Alliance pipeline.

    The curation persistent store exposes ``ontologyterm.obsolete`` as a
    boolean column, so the value reaching
    ``add_node_to_ontobio_ontology_if_not_exists`` is a Python ``bool``. The
    ``deprecated`` meta flag must be set for obsolete terms, since downstream
    filtering in ``DataManager.get_annotations_for_gene`` relies on it to skip
    obsolete GO/DO terms (whose labels are literally prefixed with
    ``"obsolete "``).
    """

    def test_boolean_true_sets_deprecated_flag(self):
        ontology = Ontology()
        AllianceDataManager.add_node_to_ontobio_ontology_if_not_exists(
            term_id="GO:0032527",
            term_label="obsolete protein exit from endoplasmic reticulum",
            term_type="biological_process",
            is_obsolete=True,
            ontology=ontology,
            check_exists=False)
        meta = ontology.node("GO:0032527")["meta"]
        self.assertTrue(meta.get("deprecated", False))

    def test_boolean_false_does_not_set_deprecated_flag(self):
        ontology = Ontology()
        AllianceDataManager.add_node_to_ontobio_ontology_if_not_exists(
            term_id="GO:0000001",
            term_label="mitochondrion inheritance",
            term_type="biological_process",
            is_obsolete=False,
            ontology=ontology,
            check_exists=False)
        meta = ontology.node("GO:0000001")["meta"]
        self.assertNotIn("deprecated", meta)

    def test_string_true_sets_deprecated_flag(self):
        """Backward compatibility: string 'true'/'True' should still work."""
        for value in ("true", "True"):
            ontology = Ontology()
            AllianceDataManager.add_node_to_ontobio_ontology_if_not_exists(
                term_id="GO:0045015",
                term_label="obsolete HDEL sequence binding",
                term_type="molecular_function",
                is_obsolete=value,
                ontology=ontology,
                check_exists=False)
            meta = ontology.node("GO:0045015")["meta"]
            self.assertTrue(meta.get("deprecated", False),
                            f"deprecated flag not set for is_obsolete={value!r}")


class TestAllianceDataManagerGeneSOTerms(unittest.TestCase):
    """Tests for the SO gene type filter applied when loading genes.

    The curation persistent store's ``gene`` table also holds non-gene sequence
    features that carry a GeneSymbolSlotAnnotation, so gene loading asks the
    curation API client to restrict results by SO gene type. Three of the four
    requested terms are not under SO "gene" in the ontology, so they cannot be
    reached from the gene root and have to be requested individually.
    """

    def setUp(self):
        this_dir = os.path.split(__file__)[0]
        conf_parser = GenedescConfigParser(os.path.join(this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.manager = AllianceDataManager(config=conf_parser)
        self.mock_db = MagicMock()
        self.mock_db.get_genes_raw.return_value = [
            {"gene_id": "MGI:97490", "gene_symbol": "Pax6"},
        ]
        self.manager._db = self.mock_db

    def test_load_gene_data_requests_the_expected_so_terms(self):
        """The four include-list SO terms are passed to the client, for the right taxon."""
        self.manager.load_gene_data(species_taxon="NCBITaxon:10090")
        self.mock_db.get_genes_raw.assert_called_once_with(
            taxon_curie="NCBITaxon:10090", so_terms=GENE_SO_TERMS)

    def test_gene_so_terms_include_mgi_and_pseudogene_roots(self):
        """gene, pseudogene, gene_segment and heritable_phenotypic_marker are all requested."""
        self.assertEqual(
            GENE_SO_TERMS,
            ["SO:0000704", "SO:0000336", "SO:3000000", "SO:0001500"])

    def test_load_gene_data_maps_rows_to_gene_objects(self):
        """Rows returned by the client become Gene objects keyed by gene id."""
        self.manager.load_gene_data(species_taxon="NCBITaxon:10090")
        self.assertIn("MGI:97490", self.manager.gene_data)
        self.assertEqual(self.manager.gene_data["MGI:97490"].name, "Pax6")


if __name__ == "__main__":
    unittest.main()
