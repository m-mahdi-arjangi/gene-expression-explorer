import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from build_copper_homeostasis_tables import (  # noqa: E402
    ALL_SAMPLES,
    CANDIDATES,
    CONTROL_SAMPLES,
    DROUGHT_SAMPLES,
    REGULON_COLUMNS,
    SAMPLE_LEVEL_COLUMNS,
    build_copper_homeostasis_tables,
    tables_match,
)


class BuildCopperHomeostasisTablesTests(unittest.TestCase):
    def setUp(self):
        self.limma = pd.DataFrame(
            [
                {"gene_id": "AT1G71200", "logFC": 0.2, "adj.P.Val": 0.03, "regulation": "upregulated"},
                {"gene_id": "AT3G46900", "logFC": 0.8, "adj.P.Val": 0.01, "regulation": "upregulated"},
                {"gene_id": "AT1G08830", "logFC": -1.2, "adj.P.Val": 0.02, "regulation": "downregulated"},
            ]
        )
        expression_rows = []
        for index, (gene_id, _symbol) in enumerate(CANDIDATES):
            values = [10 + index, 20 + index, 30 + index, 20 + index, 30 + index, 40 + index]
            expression_rows.append({"Geneid": f"{gene_id}.1", **dict(zip(ALL_SAMPLES, values))})
        self.expression = pd.DataFrame(expression_rows)
        self.table_s5 = pd.DataFrame(
            [
                {"Gene": "AT1G71200", "baseMean": 9.5, "Log2-FC": 4.0, "P (adj.)": 0.004},
                {"Gene": "AT3G46900", "baseMean": 48.2, "Log2-FC": 1.7, "P (adj.)": 0.0001},
                {"Gene": "AT1G08830", "baseMean": 2682.0, "Log2-FC": -1.4, "P (adj.)": 1e-17},
            ]
        )

    def build(self):
        return build_copper_homeostasis_tables(self.limma, self.expression, self.table_s5)

    def test_returns_fifteen_rows_in_both_tables(self):
        regulon, sample_level = self.build()
        self.assertEqual(len(regulon), 15)
        self.assertEqual(len(sample_level), 15)

    def test_preserves_curated_candidate_order(self):
        regulon, sample_level = self.build()
        expected_ids = [gene_id for gene_id, _ in CANDIDATES]
        self.assertEqual(regulon["gene_id"].tolist(), expected_ids)
        self.assertEqual(sample_level["gene_id"].tolist(), expected_ids)

    def test_returns_expected_column_order(self):
        regulon, sample_level = self.build()
        self.assertEqual(tuple(regulon.columns), REGULON_COLUMNS)
        self.assertEqual(tuple(sample_level.columns), SAMPLE_LEVEL_COLUMNS)

    def test_uses_control_and_drought_sample_groups(self):
        regulon, _ = self.build()
        first = regulon.iloc[0]
        self.assertAlmostEqual(first["GSE273926_control_mean"], np.mean([10, 20, 30]))
        self.assertAlmostEqual(first["GSE273926_drought_mean"], np.mean([20, 30, 40]))
        self.assertAlmostEqual(first["GSE273926_descriptive_ratio"], 1.5)

    def test_extracts_limma_statistics_for_candidates_that_are_present(self):
        regulon, _ = self.build()
        row = regulon.loc[regulon["gene_id"] == "AT1G71200"].iloc[0]
        self.assertAlmostEqual(row["GSE24177_logFC"], 0.2)
        self.assertAlmostEqual(row["GSE24177_FDR"], 0.03)
        self.assertEqual(row["GSE24177_regulation"], "upregulated")

    def test_marks_candidates_absent_from_limma_as_not_found(self):
        regulon, _ = self.build()
        row = regulon.loc[regulon["gene_id"] == "AT5G18830"].iloc[0]
        self.assertEqual(row["GSE24177_regulation"], "not_found")
        self.assertTrue(pd.isna(row["GSE24177_logFC"]))
        self.assertTrue(pd.isna(row["GSE24177_FDR"]))

    def test_extracts_published_table_s5_statistics(self):
        regulon, sample_level = self.build()
        comparison = regulon.loc[regulon["gene_id"] == "AT1G71200"].iloc[0]
        samples = sample_level.loc[sample_level["gene_id"] == "AT1G71200"].iloc[0]
        self.assertTrue(comparison["GSE273926_in_published_DEG_list"])
        self.assertAlmostEqual(comparison["GSE273926_published_log2FC"], 4.0)
        self.assertAlmostEqual(samples["article_baseMean"], 9.5)
        self.assertAlmostEqual(samples["article_log2FC"], 4.0)
        self.assertAlmostEqual(samples["article_FDR"], 0.004)

    def test_marks_candidates_absent_from_published_table_as_not_listed(self):
        regulon, sample_level = self.build()
        comparison = regulon.loc[regulon["gene_id"] == "AT5G18830"].iloc[0]
        samples = sample_level.loc[sample_level["gene_id"] == "AT5G18830"].iloc[0]
        self.assertFalse(comparison["GSE273926_in_published_DEG_list"])
        self.assertTrue(pd.isna(comparison["GSE273926_published_log2FC"]))
        self.assertTrue(pd.isna(samples["article_baseMean"]))

    def test_sample_level_table_preserves_six_sample_values(self):
        _, sample_level = self.build()
        row = sample_level.iloc[0]
        self.assertEqual([row[sample] for sample in CONTROL_SAMPLES], [10, 20, 30])
        self.assertEqual([row[sample] for sample in DROUGHT_SAMPLES], [20, 30, 40])

    def test_does_not_mutate_input_dataframes(self):
        limma_before = self.limma.copy(deep=True)
        expression_before = self.expression.copy(deep=True)
        s5_before = self.table_s5.copy(deep=True)
        self.build()
        pd.testing.assert_frame_equal(self.limma, limma_before)
        pd.testing.assert_frame_equal(self.expression, expression_before)
        pd.testing.assert_frame_equal(self.table_s5, s5_before)

    def test_raises_when_a_required_expression_sample_is_missing(self):
        expression = self.expression.drop(columns=[ALL_SAMPLES[0]])
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            build_copper_homeostasis_tables(self.limma, expression, self.table_s5)

    def test_raises_when_a_candidate_is_missing_from_expression_matrix(self):
        expression = self.expression.iloc[1:].copy()
        with self.assertRaisesRegex(ValueError, "missing candidate genes"):
            build_copper_homeostasis_tables(self.limma, expression, self.table_s5)

    def test_raises_when_expression_candidate_is_duplicated(self):
        expression = pd.concat([self.expression, self.expression.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "duplicate candidate gene IDs"):
            build_copper_homeostasis_tables(self.limma, expression, self.table_s5)

    def test_raises_when_expression_sample_value_is_non_numeric(self):
        expression = self.expression.copy()
        expression[ALL_SAMPLES[0]] = expression[ALL_SAMPLES[0]].astype(object)
        expression.loc[0, ALL_SAMPLES[0]] = "not-a-number"
        with self.assertRaisesRegex(ValueError, "non-numeric"):
            build_copper_homeostasis_tables(self.limma, expression, self.table_s5)

    def test_raises_when_expression_sample_value_is_negative(self):
        expression = self.expression.copy()
        expression.loc[0, ALL_SAMPLES[0]] = -1
        with self.assertRaisesRegex(ValueError, "negative"):
            build_copper_homeostasis_tables(self.limma, expression, self.table_s5)

    def test_raises_when_candidate_is_duplicated_in_limma(self):
        limma = pd.concat([self.limma, self.limma.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "Limma results contain duplicate"):
            build_copper_homeostasis_tables(limma, self.expression, self.table_s5)

    def test_raises_when_candidate_is_duplicated_in_table_s5(self):
        s5 = pd.concat([self.table_s5, self.table_s5.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "Published Table S5 U contains duplicate"):
            build_copper_homeostasis_tables(self.limma, self.expression, s5)

    def test_raises_when_limma_required_column_is_missing(self):
        limma = self.limma.drop(columns=["adj.P.Val"])
        with self.assertRaisesRegex(ValueError, "Limma results is missing required columns"):
            build_copper_homeostasis_tables(limma, self.expression, self.table_s5)

    def test_raises_when_table_s5_required_column_is_missing(self):
        s5 = self.table_s5.drop(columns=["P (adj.)"])
        with self.assertRaisesRegex(ValueError, "Published Table S5 U is missing required columns"):
            build_copper_homeostasis_tables(self.limma, self.expression, s5)

    def test_tables_match_identical_tables(self):
        regulon, _ = self.build()
        observed = regulon.copy()
        self.assertTrue(
            tables_match(
                regulon,
                observed,
                (
                    "GSE24177_logFC",
                    "GSE24177_FDR",
                    "GSE273926_control_mean",
                    "GSE273926_drought_mean",
                    "GSE273926_descriptive_ratio",
                    "GSE273926_published_log2FC",
                    "GSE273926_published_FDR",
                ),
            )
        )

    def test_tables_do_not_match_when_numeric_value_changes(self):
        regulon, _ = self.build()
        observed = regulon.copy()
        observed.loc[0, "GSE273926_control_mean"] += 1
        self.assertFalse(
            tables_match(
                regulon,
                observed,
                (
                    "GSE273926_control_mean",
                    "GSE273926_drought_mean",
                    "GSE273926_descriptive_ratio",
                ),
            )
        )

    def test_tables_do_not_match_when_column_order_changes(self):
        regulon, _ = self.build()
        observed = regulon[ list(reversed(regulon.columns)) ]
        self.assertFalse(tables_match(regulon, observed, ()))

    def test_tables_do_not_match_when_row_order_changes(self):
        regulon, _ = self.build()
        observed = regulon.iloc[::-1].reset_index(drop=True)
        self.assertFalse(
            tables_match(
                regulon,
                observed,
                REGULON_NUMERIC_COLUMNS_FOR_TEST,
            )
        )


REGULON_NUMERIC_COLUMNS_FOR_TEST = (
    "GSE24177_logFC",
    "GSE24177_FDR",
    "GSE273926_control_mean",
    "GSE273926_drought_mean",
    "GSE273926_descriptive_ratio",
    "GSE273926_published_log2FC",
    "GSE273926_published_FDR",
)


if __name__ == "__main__":
    unittest.main()
