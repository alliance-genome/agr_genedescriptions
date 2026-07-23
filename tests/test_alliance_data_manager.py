import logging
import os
import sys
import unittest

from ontobio import Ontology

sys.path.insert(0, os.path.join(os.path.split(__file__)[0], os.path.pardir, "pipelines", "alliance"))

from alliance_data_manager import AllianceDataManager  # noqa: E402

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


if __name__ == "__main__":
    unittest.main()
