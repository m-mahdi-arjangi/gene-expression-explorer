"""Tests for the curated ROS/metal gene-table reconstruction workflow."""

import unittest

import pandas as pd

from src.build_ros_metal_gene_table import (
    OUTPUT_COLUMNS,
    SELECTED_GENE_IDS,
    WorkflowValidationError,
    build_ros_metal_gene_table,
    tables_match,
)


class BuildRosMetalGeneTableTests(unittest.TestCase):
    def setUp(self):
        rows = []
        # Reverse order deliberately; the builder must return sorted gene IDs.
        for index, gene_id in enumerate(reversed(SELECTED_GENE_IDS), start=1):
            rows.append(
                {
                    "Gene": gene_id,
                    "gene_description": f"Description {gene_id}",
                    "logFC_GSE24177": -0.1 * index,
                    "FDR_GSE24177": 0.01 * index,
                    "logFC_GSE273926": -0.2 * index,
                    "FDR_GSE273926": 0.005 * index,
                    "regulation_GSE24177": "downregulated",
                    "regulation_GSE273926": "downregulated",
                    "direction_relation": "same_direction",
                }
            )
        rows.append(
            {
                "Gene": "AT9G99990",
                "gene_description": "Unselected gene",
                "logFC_GSE24177": 1.0,
                "FDR_GSE24177": 0.5,
                "logFC_GSE273926": 1.0,
                "FDR_GSE273926": 0.5,
                "regulation_GSE24177": "not_significant",
                "regulation_GSE273926": "not_significant",
                "direction_relation": "same_direction",
            }
        )
        self.source = pd.DataFrame(rows)

    def test_selects_exactly_the_eight_curated_genes(self):
        result = build_ros_metal_gene_table(self.source)
        self.assertEqual(len(result), 8)
        self.assertEqual(set(result["Gene"]), set(SELECTED_GENE_IDS))

    def test_returns_expected_column_order(self):
        result = build_ros_metal_gene_table(self.source)
        self.assertEqual(list(result.columns), list(OUTPUT_COLUMNS))

    def test_sorts_rows_by_gene_id(self):
        result = build_ros_metal_gene_table(self.source)
        self.assertEqual(result["Gene"].tolist(), sorted(SELECTED_GENE_IDS))

    def test_ignores_unselected_source_genes(self):
        result = build_ros_metal_gene_table(self.source)
        self.assertNotIn("AT9G99990", set(result["Gene"]))

    def test_preserves_selected_source_values(self):
        result = build_ros_metal_gene_table(self.source)
        source_row = self.source.loc[
            self.source["Gene"] == "AT1G08830"
        ].iloc[0]
        result_row = result.loc[result["Gene"] == "AT1G08830"].iloc[0]
        self.assertEqual(
            result_row["gene_description"], source_row["gene_description"]
        )
        self.assertEqual(
            result_row["logFC_GSE24177"], source_row["logFC_GSE24177"]
        )
        self.assertEqual(
            result_row["FDR_GSE273926"], source_row["FDR_GSE273926"]
        )

    def test_does_not_mutate_input_dataframe(self):
        original = self.source.copy(deep=True)
        build_ros_metal_gene_table(self.source)
        pd.testing.assert_frame_equal(self.source, original)

    def test_raises_when_a_selected_gene_is_missing(self):
        incomplete = self.source.loc[self.source["Gene"] != "AT1G08830"]
        with self.assertRaisesRegex(WorkflowValidationError, "AT1G08830"):
            build_ros_metal_gene_table(incomplete)

    def test_raises_when_a_selected_gene_is_duplicated(self):
        duplicate = pd.concat(
            [self.source, self.source.loc[self.source["Gene"] == "AT1G08830"]],
            ignore_index=True,
        )
        with self.assertRaisesRegex(WorkflowValidationError, "Duplicate"):
            build_ros_metal_gene_table(duplicate)

    def test_raises_when_a_required_column_is_missing(self):
        incomplete = self.source.drop(columns=["direction_relation"])
        with self.assertRaisesRegex(WorkflowValidationError, "direction_relation"):
            build_ros_metal_gene_table(incomplete)

    def test_raises_when_a_numeric_value_is_not_numeric(self):
        invalid = self.source.copy()
        invalid["logFC_GSE24177"] = invalid["logFC_GSE24177"].astype(object)
        invalid.loc[invalid["Gene"] == "AT1G08830", "logFC_GSE24177"] = "bad"
        with self.assertRaisesRegex(WorkflowValidationError, "logFC_GSE24177"):
            build_ros_metal_gene_table(invalid)

    def test_raises_when_a_numeric_value_is_missing(self):
        invalid = self.source.copy()
        invalid.loc[invalid["Gene"] == "AT1G08830", "FDR_GSE24177"] = None
        with self.assertRaisesRegex(WorkflowValidationError, "FDR_GSE24177"):
            build_ros_metal_gene_table(invalid)

    def test_raises_when_fdr_is_outside_zero_to_one(self):
        invalid = self.source.copy()
        invalid.loc[invalid["Gene"] == "AT1G08830", "FDR_GSE24177"] = 1.2
        with self.assertRaisesRegex(WorkflowValidationError, r"outside \[0, 1\]"):
            build_ros_metal_gene_table(invalid)

    def test_tables_match_after_row_order_normalization(self):
        expected = build_ros_metal_gene_table(self.source)
        observed = expected.sample(frac=1, random_state=7).reset_index(drop=True)
        self.assertTrue(tables_match(expected, observed))

    def test_tables_do_not_match_when_a_value_changes(self):
        expected = build_ros_metal_gene_table(self.source)
        observed = expected.copy()
        observed.loc[0, "logFC_GSE24177"] += 0.01
        self.assertFalse(tables_match(expected, observed))

    def test_tables_do_not_match_when_columns_differ(self):
        expected = build_ros_metal_gene_table(self.source)
        observed = expected.drop(columns=["direction_relation"])
        self.assertFalse(tables_match(expected, observed))


if __name__ == "__main__":
    unittest.main()
