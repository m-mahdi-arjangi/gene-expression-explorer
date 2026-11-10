"""Unit tests for reproducible cross-study drought analysis."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from math import isclose

import numpy as np
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

import cross_study_analysis as analysis


class CrossStudyAnalysisTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        base = Path(self.temp_dir.name)

        self.expression_path = base / "expression.csv"
        self.limma_path = base / "limma.csv"
        self.table_s5_path = base / "table_s5.csv"
        self.output_dir = base / "outputs"

        self.default_expression_ids = [
            "G1", "G2", "G3", "G4", "G5",
            "G6", "G7", "G8", "G9", "G10",
            "ATCG00020",
            "EXPR_ONLY",
        ]

        self.default_limma = pd.DataFrame(
            [
                ("G1", 1.2, 0.01),
                ("G2", -1.1, 0.02),
                ("G3", 0.8, 0.03),
                ("G4", 0.0, 0.04),
                ("G5", -0.9, 0.01),
                ("G6", 0.2, 0.20),
                ("G7", -0.3, 0.30),
                ("G8", 0.4, 0.50),
                ("ATCG00020", 0.1, 0.60),
                ("LIMMA_ONLY", 0.2, 0.70),
            ],
            columns=["gene_id", "logFC", "adj.P.Val"],
        )

        self.default_limma["probe_id"] = (
            self.default_limma["gene_id"] + "_at"
        )
        self.default_limma["AveExpr"] = 8.0
        self.default_limma["t"] = 2.0
        self.default_limma["P.Value"] = 0.01
        self.default_limma["B"] = 1.0
        self.default_limma["locus_tag"] = (
            self.default_limma["gene_id"]
        )
        self.default_limma["gene_description"] = (
            "Synthetic test gene"
        )

        self.default_limma = self.default_limma[
            [
                "probe_id",
                "logFC",
                "AveExpr",
                "t",
                "P.Value",
                "adj.P.Val",
                "B",
                "locus_tag",
                "gene_id",
                "gene_description",
            ]
        ]

        self.default_s5 = pd.DataFrame(
            [
                ("G1", 0.6, 0.001),
                ("G2", -0.7, 0.002),
                ("G3", -0.6, 0.003),
                ("G4", 0.5, 0.004),
                ("G5", 0.8, 0.005),
                ("G6", -0.4, 0.010),
                ("EXPR_ONLY", 0.4, 0.020),
                ("S5_ONLY", 0.3, 0.030),
            ],
            columns=["Gene", "Log2-FC", "P (adj.)"],
        )

        patches = [
            patch.object(
                analysis, "EXPRESSION_FILE", self.expression_path
            ),
            patch.object(
                analysis, "LIMMA_FILE", self.limma_path
            ),
            patch.object(
                analysis, "TABLE_S5_FILE", self.table_s5_path
            ),
            patch.object(
                analysis, "OUTPUT_DIR", self.output_dir
            ),
        ]

        for item in patches:
            item.start()
            self.addCleanup(item.stop)

        self.write_inputs()

    def write_inputs(
        self,
        expression_ids=None,
        limma=None,
        table_s5=None,
    ):
        """Create isolated CSV inputs for each test."""

        if expression_ids is None:
            expression_ids = self.default_expression_ids

        if limma is None:
            limma = self.default_limma.copy()

        if table_s5 is None:
            table_s5 = self.default_s5.copy()

        pd.DataFrame(
            {"gene_id": expression_ids}
        ).to_csv(self.expression_path, index=False)

        limma.to_csv(self.limma_path, index=False)
        table_s5.to_csv(self.table_s5_path, index=False)

    # --------------------------------------------------
    # Identifier and input validation
    # --------------------------------------------------

    def test_normalize_ids_removes_blank_and_invalid_values(self):
        values = pd.Series(
            [
                " G1 ",
                None,
                "",
                "nan",
                "NONE",
                "<NA>",
                "ATCG00020",
            ]
        )

        result = analysis.normalize_ids(values).tolist()

        self.assertEqual(
            result,
            ["G1", "ATCG00020"],
        )

    def test_require_columns_rejects_missing_columns(self):
        dataframe = pd.DataFrame(
            {"gene_id": ["G1"]}
        )

        with self.assertRaisesRegex(
            ValueError,
            "missing columns",
        ):
            analysis.require_columns(
                dataframe,
                {"gene_id", "logFC"},
                "Test input",
            )

    def test_rejects_duplicate_expression_ids(self):
        duplicated_ids = (
            self.default_expression_ids + ["G1"]
        )

        self.write_inputs(
            expression_ids=duplicated_ids
        )

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate gene IDs",
        ):
            analysis.calculate_results()

    def test_rejects_duplicate_significant_limma_ids(self):
        limma = self.default_limma.copy()

        duplicate = limma.loc[
            limma["gene_id"] == "G1"
        ].copy()

        limma = pd.concat(
            [limma, duplicate],
            ignore_index=True,
        )

        self.write_inputs(limma=limma)

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate significant gene IDs",
        ):
            analysis.calculate_results()

    def test_rejects_duplicate_table_s5_ids(self):
        table_s5 = self.default_s5.copy()

        duplicate = table_s5.loc[
            table_s5["Gene"] == "G1"
        ].copy()

        table_s5 = pd.concat(
            [table_s5, duplicate],
            ignore_index=True,
        )

        self.write_inputs(table_s5=table_s5)

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate significant gene IDs",
        ):
            analysis.calculate_results()

    def test_rejects_missing_limma_file(self):
        self.limma_path.unlink()

        with self.assertRaises(FileNotFoundError):
            analysis.calculate_results()

    def test_rejects_missing_table_s5_column(self):
        table_s5 = self.default_s5.drop(
            columns=["Log2-FC"]
        )

        self.write_inputs(table_s5=table_s5)

        with self.assertRaisesRegex(
            ValueError,
            "missing columns",
        ):
            analysis.calculate_results()

    def test_rejects_nonfinite_limma_logfc(self):
        limma = self.default_limma.copy()
        limma.loc[
            limma["gene_id"] == "G1", "logFC"
        ] = np.inf

        self.write_inputs(limma=limma)

        with self.assertRaisesRegex(
            ValueError,
            "Non-finite Limma",
        ):
            analysis.calculate_results()

    def test_rejects_nonfinite_table_s5_logfc(self):
        table_s5 = self.default_s5.copy()
        table_s5.loc[
            table_s5["Gene"] == "G1", "Log2-FC"
        ] = np.inf

        self.write_inputs(table_s5=table_s5)

        with self.assertRaisesRegex(
            ValueError,
            "Non-finite Table S5",
        ):
            analysis.calculate_results()

    # --------------------------------------------------
    # Common universe and overlap calculations
    # --------------------------------------------------

    def test_calculates_expected_universe_and_overlap(self):
        outputs, metrics = analysis.calculate_results()

        self.assertEqual(metrics["universe"], 9)
        self.assertEqual(metrics["limma_significant"], 5)
        self.assertEqual(
            metrics["table_s5_significant_total"], 8
        )
        self.assertEqual(
            metrics["table_s5_inside_universe"], 6
        )
        self.assertEqual(
            metrics["table_s5_outside_universe"], 2
        )
        self.assertEqual(metrics["shared"], 5)

        universe = outputs[
            "common_tested_gene_universe.csv"
        ][0]

        self.assertEqual(len(universe), 9)
        self.assertIn(
            "ATCG00020",
            set(universe["gene_id"]),
        )
        self.assertNotIn(
            "EXPR_ONLY",
            set(universe["gene_id"]),
        )
        self.assertNotIn(
            "LIMMA_ONLY",
            set(universe["gene_id"]),
        )

    def test_calculates_expected_direction_groups(self):
        outputs, metrics = analysis.calculate_results()

        self.assertEqual(metrics["same_direction"], 2)
        self.assertEqual(metrics["opposite_direction"], 2)
        self.assertEqual(metrics["unresolved"], 1)

        conserved = outputs[
            "conserved_drought_signature_45_genes.csv"
        ][0]

        upregulated = outputs[
            "shared_upregulated_drought_genes_29.csv"
        ][0]

        downregulated = outputs[
            "shared_downregulated_drought_genes_16.csv"
        ][0]

        discordant = outputs[
            "discordant_drought_genes_7.csv"
        ][0]

        self.assertEqual(
            set(conserved["Gene"]),
            {"G1", "G2"},
        )
        self.assertEqual(
            set(upregulated["Gene"]),
            {"G1"},
        )
        self.assertEqual(
            set(downregulated["Gene"]),
            {"G2"},
        )
        self.assertEqual(
            set(discordant["Gene"]),
            {"G3", "G5"},
        )

        self.assertEqual(
            len(conserved)
            + len(discordant)
            + metrics["unresolved"],
            metrics["shared"],
        )

    def test_zero_logfc_is_classified_as_unresolved(self):
        outputs, _ = analysis.calculate_results()

        shared = outputs[
            "significant_drought_genes_shared.csv"
        ][0]

        gene = shared.loc[
            shared["Gene"] == "G4"
        ].iloc[0]

        self.assertEqual(
            gene["regulation_GSE24177"],
            "unresolved",
        )
        self.assertEqual(
            gene["direction_relation"],
            "unresolved",
        )

    # --------------------------------------------------
    # Statistical calculations
    # --------------------------------------------------

    def test_calculates_expected_enrichment_statistics(self):
        outputs, _ = analysis.calculate_results()

        stats = outputs[
            "cross_study_statistical_validation.csv"
        ][0].iloc[0]

        self.assertEqual(
            stats["common_tested_universe"], 9
        )
        self.assertEqual(
            stats["significant_GSE24177"], 5
        )
        self.assertEqual(
            stats["significant_GSE273926"], 6
        )
        self.assertEqual(
            stats["observed_overlap"], 5
        )

        self.assertTrue(
            isclose(
                stats["expected_overlap"],
                10 / 3,
                rel_tol=1e-12,
            )
        )

        self.assertTrue(
            isclose(
                stats["fold_enrichment"],
                1.5,
                rel_tol=1e-12,
            )
        )

        self.assertTrue(
            isclose(
                stats["hypergeometric_p"],
                1 / 21,
                rel_tol=1e-10,
            )
        )

    def test_calculates_direction_contingency_table(self):
        outputs, _ = analysis.calculate_results()

        stats = outputs[
            "cross_study_statistical_validation.csv"
        ][0].iloc[0]

        self.assertEqual(stats["up_up"], 1)
        self.assertEqual(stats["up_down"], 1)
        self.assertEqual(stats["down_up"], 1)
        self.assertEqual(stats["down_down"], 1)

        self.assertTrue(
            isclose(
                stats["direction_fisher_odds_ratio"],
                1.0,
            )
        )
        self.assertTrue(
            isclose(
                stats["direction_fisher_p"],
                1.0,
            )
        )

    def test_overlap_summary_counts_are_consistent(self):
        outputs, _ = analysis.calculate_results()

        summary = outputs[
            "cross_study_overlap_summary.csv"
        ][0].iloc[0]

        self.assertEqual(
            summary["shared_significant_genes"], 5
        )
        self.assertEqual(summary["same_direction"], 2)
        self.assertEqual(summary["opposite_direction"], 2)
        self.assertEqual(summary["unresolved"], 1)

        self.assertTrue(
            isclose(
                summary["same_direction_pct_of_overlap"],
                40.0,
            )
        )

    # --------------------------------------------------
    # Comparison with existing output files
    # --------------------------------------------------

    def test_existing_outputs_match_rebuilt_results(self):
        outputs, _ = analysis.calculate_results()

        self.output_dir.mkdir(parents=True)

        for filename, (dataframe, _) in outputs.items():
            dataframe.to_csv(
                self.output_dir / filename,
                index=False,
            )

        mismatches = analysis.compare_existing_outputs(
            outputs
        )

        self.assertEqual(mismatches, [])

    def test_comparison_detects_changed_output_value(self):
        outputs, _ = analysis.calculate_results()

        self.output_dir.mkdir(parents=True)

        for filename, (dataframe, _) in outputs.items():
            dataframe.to_csv(
                self.output_dir / filename,
                index=False,
            )

        summary_path = (
            self.output_dir
            / "cross_study_overlap_summary.csv"
        )

        saved_summary = pd.read_csv(summary_path)
        saved_summary.loc[0, "same_direction"] = 999

        saved_summary.to_csv(
            summary_path,
            index=False,
        )

        mismatches = analysis.compare_existing_outputs(
            outputs
        )

        self.assertIn(
            "cross_study_overlap_summary.csv",
            mismatches,
        )


if __name__ == "__main__":
    unittest.main()
