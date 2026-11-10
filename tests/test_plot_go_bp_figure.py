import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

import plot_go_bp_figure as plot


class PlotFigureValidationTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.base = Path(self.temp_dir.name)

        self.membership_path = self.base / "GO_BP_term_gene_membership.csv"

        self.membership = pd.DataFrame([
            {
                "gene_set": "Conserved 45",
                "GO_ID": "GO:0000001",
                "GO_term": "test process one",
                "adjusted_p_value": 0.01,
                "gene_id": "AT1G00010",
                "evidence_codes": "IEA",
            },
            {
                "gene_set": "Conserved 45",
                "GO_ID": "GO:0000001",
                "GO_term": "test process one",
                "adjusted_p_value": 0.01,
                "gene_id": "AT2G00020",
                "evidence_codes": "IMP",
            },
            {
                "gene_set": "Conserved 45",
                "GO_ID": "GO:0000002",
                "GO_term": "test process two",
                "adjusted_p_value": 0.02,
                "gene_id": "AT2G00020",
                "evidence_codes": "IEA",
            },
        ])

        self.enrichment = pd.DataFrame([
            {
                "native": "GO:0000001",
                "name": "test process one",
                "p_value": 0.01,
                "intersection_size": 2,
                "significant": True,
                "source": "GO:BP",
            },
            {
                "native": "GO:0000002",
                "name": "test process two",
                "p_value": 0.02,
                "intersection_size": 1,
                "significant": True,
                "source": "GO:BP",
            },
        ])

        self.enrichment_path = (
            self.base / "conserved_45_GO_BP_enrichment.csv"
        )

        self.membership.to_csv(
            self.membership_path,
            index=False,
        )

        self.enrichment.to_csv(
            self.enrichment_path,
            index=False,
        )

        self.data_patch = patch.object(
            plot,
            "DATA_DIR",
            self.base,
        )

        self.membership_patch = patch.object(
            plot,
            "MEMBERSHIP_FILE",
            self.membership_path,
        )

        self.data_patch.start()
        self.membership_patch.start()

    def tearDown(self):
        self.data_patch.stop()
        self.membership_patch.stop()
        self.temp_dir.cleanup()

    def test_loads_valid_membership_table(self):
        result = plot.load_membership()

        self.assertEqual(len(result), 3)
        self.assertEqual(
            result["gene_id"].nunique(),
            2,
        )

    def test_rejects_missing_membership_column(self):
        self.membership.drop(
            columns=["evidence_codes"]
        ).to_csv(
            self.membership_path,
            index=False,
        )

        with self.assertRaisesRegex(
            ValueError,
            "missing columns",
        ):
            plot.load_membership()

    def test_rejects_duplicate_gene_term_memberships(self):
        duplicated = pd.concat(
            [self.membership, self.membership.iloc[[0]]],
            ignore_index=True,
        )

        duplicated.to_csv(
            self.membership_path,
            index=False,
        )

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate gene-term memberships",
        ):
            plot.load_membership()

    def test_loads_and_validates_enrichment_results(self):
        membership = plot.load_membership()

        result = plot.load_enrichment(
            "Conserved 45",
            plot.GENE_SETS["Conserved 45"],
            membership,
        )

        self.assertEqual(len(result), 2)
        self.assertEqual(
            result.iloc[0]["native"],
            "GO:0000001",
        )

    def test_rejects_missing_membership_term(self):
        membership = self.membership.loc[
            self.membership["GO_ID"] == "GO:0000001"
        ].copy()

        with self.assertRaisesRegex(
            ValueError,
            "GO-term mismatch",
        ):
            plot.load_enrichment(
                "Conserved 45",
                plot.GENE_SETS["Conserved 45"],
                membership,
            )

    def test_rejects_intersection_size_mismatch(self):
        membership = self.membership.loc[
            ~(
                (self.membership["GO_ID"] == "GO:0000001")
                & (self.membership["gene_id"] == "AT2G00020")
            )
        ].copy()

        with self.assertRaisesRegex(
            ValueError,
            "Intersection-size mismatch",
        ):
            plot.load_enrichment(
                "Conserved 45",
                plot.GENE_SETS["Conserved 45"],
                membership,
            )

    def test_rejects_adjusted_p_value_mismatch(self):
        membership = self.membership.copy()

        membership.loc[
            membership["GO_ID"] == "GO:0000001",
            "adjusted_p_value",
        ] = 0.03

        with self.assertRaisesRegex(
            ValueError,
            "Adjusted p-value mismatch",
        ):
            plot.load_enrichment(
                "Conserved 45",
                plot.GENE_SETS["Conserved 45"],
                membership,
            )

    def test_rejects_non_significant_only_enrichment(self):
        enrichment = self.enrichment.copy()
        enrichment["significant"] = False

        enrichment.to_csv(
            self.enrichment_path,
            index=False,
        )

        membership = plot.load_membership()

        with self.assertRaisesRegex(
            ValueError,
            "No significant GO:BP terms",
        ):
            plot.load_enrichment(
                "Conserved 45",
                plot.GENE_SETS["Conserved 45"],
                membership,
            )


if __name__ == "__main__":
    unittest.main()
