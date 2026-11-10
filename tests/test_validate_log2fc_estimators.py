"""Unit tests for GSE273926 log2FC reconstruction and validation."""

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

from validate_log2fc_estimators import (
    DIRECTION_OUTPUT_COLUMNS,
    SAMPLE_COLUMNS,
    build_direction_summary,
    build_estimator_direction_differences,
    build_outputs,
    calculate_estimators,
    calculate_metrics,
    classify_direction,
    classify_relation,
    compare_existing_outputs,
    dataframes_equal_by_gene,
    matches_expected_scientific_correction,
    validate_table_s4,
    validate_table_s5,
)


def make_s4():
    """Create a small synthetic Table S4 with four genes."""

    genes = ["ATG00001", "ATG00002", "ATG00003", "ATG00004"]

    values = {
        "U": {
            "ww": [
                [10, 10, 10],
                [20, 20, 20],
                [10, 10, 10],
                [30, 30, 30],
            ],
            "md": [
                [20, 20, 20],
                [40, 40, 40],
                [5, 10, 15],
                [15, 15, 15],
            ],
        },
        "D": {
            "ww": [
                [10, 10],
                [20, 20],
                [10, 10],
                [20, 20],
            ],
            "md": [
                [30, 30],
                [40, 40],
                [5, 15],
                [40, 40],
            ],
        },
        "FD": {
            "ww": [
                [10, 10],
                [20, 20],
                [10, 10],
                [20, 20],
            ],
            "md": [
                [40, 40],
                [10, 30],
                [5, 15],
                [10, 10],
            ],
        },
    }

    data = {"Gene": genes}

    for group, conditions in values.items():
        for condition in ("ww", "md"):
            columns = SAMPLE_COLUMNS[group][condition]
            rows = conditions[condition]

            for sample_index, column in enumerate(columns):
                data[column] = [
                    row[sample_index]
                    for row in rows
                ]

    return pd.DataFrame(data)


def make_s5(group):
    """Create a small synthetic Table S5 for one sample type."""

    genes = [
        "ATG00001",
        "ATG00002",
        "ATG00003",
        "ATG00004",
    ]

    logfc = {
        "U": [1.2, 0.8, -0.6, -1.1],
        "D": [1.4, 1.0, -0.5, 0.2],
        "FD": [1.6, -0.5, -0.8, 0.1],
    }[group]

    data = {
        "Gene": genes,
        "baseMean": [100.0, 200.0, 300.0, 400.0],
        "Log2-FC": logfc,
        "lfcSE": [0.2, 0.3, 0.25, 0.4],
        "Stat.": [3.0, 2.5, -2.4, 2.0],
        "P-value": [0.001, 0.002, 0.003, 0.004],
        "P (adj.)": [0.01, 0.02, 0.03, 0.04],
    }

    if group in ("D", "FD"):
        data["DE in U samples"] = ["Yes", "Yes", "Yes", "No"]

    return pd.DataFrame(data)


class SyntheticOutputTests(unittest.TestCase):

    def setUp(self):
        self.s4 = make_s4()
        self.u = make_s5("U")
        self.d = make_s5("D")
        self.fd = make_s5("FD")

    def test_build_outputs_returns_all_expected_files(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        self.assertEqual(len(outputs), 12)
        self.assertIn("log2FC_estimator_validation.csv", outputs)
        self.assertIn("U_701_direction_from_tableS4.csv", outputs)
        self.assertIn("FD_estimator_direction_differences.csv", outputs)

    def test_build_outputs_preserves_all_u_genes(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        direction = outputs["U_701_direction_from_tableS4.csv"]

        self.assertEqual(len(direction), 4)
        self.assertEqual(
            direction["Gene"].tolist(),
            self.u["Gene"].tolist(),
        )

    def test_build_outputs_filters_d_and_fd_by_de_in_u(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        self.assertEqual(
            len(outputs["D_published_vs_estimated_log2FC.csv"]),
            3,
        )
        self.assertEqual(
            len(outputs["FD_published_vs_estimated_log2FC.csv"]),
            3,
        )

    def test_build_outputs_calculates_expected_validation_rows(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        validation = outputs["log2FC_estimator_validation.csv"]

        self.assertEqual(len(validation), 6)
        self.assertEqual(
            validation["sample_type"].tolist(),
            ["U", "U", "D", "D", "FD", "FD"],
        )
        self.assertEqual(
            validation["method"].tolist(),
            [1, 2, 1, 2, 1, 2],
        )

    def test_build_outputs_recognizes_zero_direction_as_unresolved(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        summary = outputs[
            "U_701_direction_reconstruction_summary.csv"
        ]

        fd_summary = summary.loc[
            summary["sample_type"].eq("FD")
        ].iloc[0]

        self.assertEqual(fd_summary["unresolved"], 2)

    def test_build_outputs_preserves_published_u_d_opposite_genes(self):
        outputs = build_outputs(
            self.s4,
            self.u,
            self.d,
            self.fd,
        )

        opposite = outputs["U_vs_D_opposite_direction_genes.csv"]

        self.assertIn("direction_relation",
                      outputs["article_table_s5_U_vs_D.csv"].columns)
        self.assertTrue(
            set(opposite["Gene"]).issubset(set(self.u["Gene"]))
        )


class EstimatorTests(unittest.TestCase):

    def setUp(self):
        self.s4 = make_s4()

    def test_ratio_of_means_uses_pseudocount_one(self):
        result = calculate_estimators(self.s4)

        expected = np.log2(21 / 11)

        self.assertAlmostEqual(
            result.loc[0, "U_S4_log2FC"],
            expected,
            places=12,
        )

    def test_mean_log2_estimator_uses_pseudocount_one(self):
        result = calculate_estimators(self.s4)

        ww = np.array([10, 10, 10], dtype=float)
        md = np.array([20, 20, 20], dtype=float)

        expected = (
            np.log2(md + 1).mean()
            - np.log2(ww + 1).mean()
        )

        self.assertAlmostEqual(
            result.loc[0, "U_mean_log2FC"],
            expected,
            places=12,
        )

    def test_estimators_are_not_assumed_to_be_identical(self):
        result = calculate_estimators(self.s4)

        # Equal group means produce a zero ratio-of-means estimator,
        # while different replicate distributions can produce a
        # non-zero mean-log2 estimator.
        self.assertAlmostEqual(
            result.loc[2, "U_S4_log2FC"],
            0.0,
            places=12,
        )

        self.assertNotAlmostEqual(
            result.loc[2, "U_mean_log2FC"],
            0.0,
            places=6,
        )

    def test_estimator_output_contains_expected_columns(self):
        result = calculate_estimators(self.s4)

        self.assertEqual(
            result.columns.tolist(),
            [
                "Gene",
                "U_S4_log2FC",
                "U_mean_log2FC",
                "D_S4_log2FC",
                "D_mean_log2FC",
                "FD_S4_log2FC",
                "FD_mean_log2FC",
            ],
        )


class DirectionTests(unittest.TestCase):

    def test_classifies_positive_negative_and_zero(self):
        result = classify_direction([1.0, -1.0, 0.0])

        self.assertEqual(
            result.tolist(),
            ["up", "down", "zero"],
        )

    def test_rejects_nonfinite_direction_values(self):
        with self.assertRaises(ValueError):
            classify_direction([1.0, np.nan])

    def test_zero_direction_is_unresolved(self):
        result = classify_relation(
            [1.0, -1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0, -1.0],
        )

        self.assertEqual(
            result.tolist(),
            [
                "same_direction",
                "opposite_direction",
                "unresolved",
                "unresolved",
            ],
        )

    def test_direction_summary_counts_unresolved_separately(self):
        direction = pd.DataFrame({
            "U_article_direction": ["up", "down", "up", "down"],
            "U_S4_direction": ["up", "down", "down", "down"],
            "D_S4_direction": ["up", "zero", "up", "down"],
            "FD_S4_direction": ["up", "down", "zero", "up"],
        })

        result = build_direction_summary(direction)

        fd = result.loc[result["sample_type"].eq("FD")].iloc[0]

        self.assertEqual(fd["same_direction"], 2)
        self.assertEqual(fd["opposite_direction"], 1)
        self.assertEqual(fd["unresolved"], 1)

        d = result.loc[result["sample_type"].eq("D")].iloc[0]

        self.assertEqual(d["same_direction"], 3)
        self.assertEqual(d["opposite_direction"], 0)
        self.assertEqual(d["unresolved"], 1)

    def test_estimator_differences_include_zero_vs_nonzero(self):
        direction = pd.DataFrame({
            column: [0, 0, 0]
            for column in DIRECTION_OUTPUT_COLUMNS
        })

        direction["Gene"] = ["A", "B", "C"]
        direction["FD_S4_direction"] = ["zero", "up", "down"]
        direction["FD_mean_log_direction"] = ["down", "down", "down"]

        result = build_estimator_direction_differences(
            direction,
            "FD",
        )

        self.assertEqual(result["Gene"].tolist(), ["A", "B"])
        self.assertEqual(
            result.columns.tolist(),
            DIRECTION_OUTPUT_COLUMNS,
        )


class InputValidationTests(unittest.TestCase):

    def setUp(self):
        self.s4 = make_s4()
        self.u = make_s5("U")
        self.d = make_s5("D")

    def test_accepts_valid_table_s4(self):
        result = validate_table_s4(self.s4)
        self.assertEqual(len(result), 4)

    def test_rejects_missing_table_s4_column(self):
        invalid = self.s4.drop(columns=["U_WW_1-1"])

        with self.assertRaises(ValueError):
            validate_table_s4(invalid)

    def test_rejects_duplicate_table_s4_gene_ids(self):
        invalid = self.s4.copy()
        invalid.loc[1, "Gene"] = invalid.loc[0, "Gene"]

        with self.assertRaises(ValueError):
            validate_table_s4(invalid)

    def test_rejects_negative_table_s4_expression(self):
        invalid = self.s4.copy()
        invalid.loc[0, "U_WW_1-1"] = -1

        with self.assertRaises(ValueError):
            validate_table_s4(invalid)

    def test_rejects_nonfinite_table_s4_expression(self):
        invalid = self.s4.copy()
        invalid["U_WW_1-1"] = invalid["U_WW_1-1"].astype(float)
        invalid.loc[0, "U_WW_1-1"] = np.inf

        with self.assertRaises(ValueError):
            validate_table_s4(invalid)

    def test_accepts_valid_table_s5(self):
        result = validate_table_s5(self.u, "U")
        self.assertEqual(len(result), 4)

    def test_rejects_duplicate_table_s5_gene_ids(self):
        invalid = self.u.copy()
        invalid.loc[1, "Gene"] = invalid.loc[0, "Gene"]

        with self.assertRaises(ValueError):
            validate_table_s5(invalid, "U")

    def test_rejects_adjusted_p_value_outside_range(self):
        invalid = self.u.copy()
        invalid.loc[0, "P (adj.)"] = 1.5

        with self.assertRaises(ValueError):
            validate_table_s5(invalid, "U")

    def test_rejects_invalid_de_in_u_flag(self):
        invalid = self.d.copy()
        invalid.loc[0, "DE in U samples"] = "Maybe"

        with self.assertRaises(ValueError):
            validate_table_s5(invalid, "D")


class MetricsTests(unittest.TestCase):

    def test_calculates_metrics_and_excludes_zero_from_direction_agreement(self):
        result = calculate_metrics(
            [1.0, -1.0, 0.0],
            [2.0, 1.0, -1.0],
        )

        self.assertEqual(result["genes_compared"], 3)
        self.assertAlmostEqual(result["MAE"], 4 / 3)
        self.assertAlmostEqual(result["RMSE"], np.sqrt(2))
        self.assertAlmostEqual(result["mean_bias"], 2 / 3)
        self.assertAlmostEqual(
            result["direction_agreement_pct"],
            50.0,
        )


class DataFrameComparisonTests(unittest.TestCase):

    def test_gene_keyed_comparison_ignores_row_order(self):
        left = pd.DataFrame({
            "Gene": ["B", "A"],
            "value": [2.0, 1.0],
        })
        right = pd.DataFrame({
            "Gene": ["A", "B"],
            "value": [1.0, 2.0],
        })

        self.assertTrue(dataframes_equal_by_gene(left, right))

    def test_gene_keyed_comparison_detects_value_difference(self):
        left = pd.DataFrame({
            "Gene": ["A", "B"],
            "value": [1.0, 2.0],
        })
        right = pd.DataFrame({
            "Gene": ["A", "B"],
            "value": [1.0, 3.0],
        })

        self.assertFalse(dataframes_equal_by_gene(left, right))

    def test_gene_keyed_comparison_rejects_duplicate_ids(self):
        left = pd.DataFrame({
            "Gene": ["A", "A"],
            "value": [1.0, 2.0],
        })
        right = pd.DataFrame({
            "Gene": ["A", "B"],
            "value": [1.0, 2.0],
        })

        self.assertFalse(dataframes_equal_by_gene(left, right))

    def test_gene_keyed_comparison_rejects_different_columns(self):
        left = pd.DataFrame({
            "Gene": ["A"],
            "value": [1.0],
        })
        right = pd.DataFrame({
            "Gene": ["A"],
            "other": [1.0],
        })

        self.assertFalse(dataframes_equal_by_gene(left, right))


class ScientificCorrectionTests(unittest.TestCase):

    def test_accepts_at2g16740_direction_correction(self):
        existing = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["opposite_direction"],
        })
        rebuilt = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["unresolved"],
        })

        self.assertTrue(
            matches_expected_scientific_correction(
                "U_701_direction_from_tableS4.csv",
                existing,
                rebuilt,
            )
        )

    def test_rejects_unexpected_direction_table_difference(self):
        existing = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["opposite_direction"],
            "other": [1],
        })
        rebuilt = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["unresolved"],
            "other": [2],
        })

        self.assertFalse(
            matches_expected_scientific_correction(
                "U_701_direction_from_tableS4.csv",
                existing,
                rebuilt,
            )
        )

    def test_accepts_fd_summary_zero_correction(self):
        existing = pd.DataFrame({
            "sample_type": ["FD"],
            "genes_compared": [701],
            "same_direction": [661],
            "opposite_direction": [40],
            "unresolved": [0],
        })
        rebuilt = pd.DataFrame({
            "sample_type": ["FD"],
            "genes_compared": [701],
            "same_direction": [661],
            "opposite_direction": [39],
            "unresolved": [1],
        })

        self.assertTrue(
            matches_expected_scientific_correction(
                "U_701_direction_reconstruction_summary.csv",
                existing,
                rebuilt,
            )
        )

    def test_accepts_fd_estimator_difference_correction(self):
        existing = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["opposite_direction"],
        })
        rebuilt = pd.DataFrame({
            "Gene": ["AT2G16740"],
            "FD_direction_relation": ["unresolved"],
        })

        self.assertTrue(
            matches_expected_scientific_correction(
                "FD_estimator_direction_differences.csv",
                existing,
                rebuilt,
            )
        )


class ExistingOutputComparisonTests(unittest.TestCase):

    def setUp(self):
        self.outputs = build_outputs(
            make_s4(),
            make_s5("U"),
            make_s5("D"),
            make_s5("FD"),
        )
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

        for filename, dataframe in self.outputs.items():
            dataframe.to_csv(
                self.output_dir / filename,
                index=False,
            )

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_accepts_matching_existing_outputs(self):
        needs_update, conflicts = compare_existing_outputs(
            self.outputs,
            self.output_dir,
        )

        self.assertEqual(needs_update, [])
        self.assertEqual(conflicts, [])

    def test_detects_unexpected_existing_output_difference(self):
        filename = "U_published_vs_estimated_log2FC.csv"
        path = self.output_dir / filename

        existing = pd.read_csv(path)
        existing.loc[0, "Log2-FC"] += 0.5
        existing.to_csv(path, index=False)

        needs_update, conflicts = compare_existing_outputs(
            self.outputs,
            self.output_dir,
        )

        self.assertIn(filename, conflicts)
        self.assertNotIn(filename, needs_update)


if __name__ == "__main__":
    unittest.main()
