import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

from extract_go_term_membership import (
    GENE_SETS,
    OUTPUT_COLUMNS,
    build_membership_table,
    compare_membership_tables,
    extract_membership_rows,
)


def make_response():
    """Create a small, valid g:Profiler-like response for testing."""

    genes = ["AT1G00010", "AT2G00020", "AT3G00030"]

    return {
        "meta": {
            "query_metadata": {
                "organism": "athaliana",
                "significance_threshold_method": "fdr",
                "user_threshold": 0.05,
                "sources": ["GO:BP"],
                "queries": {
                    "query_1": genes,
                },
            },
            "genes_metadata": {
                "query": {
                    "query_1": {
                        "ensgs": genes,
                    },
                },
            },
        },
        "result": [
            {
                "source": "GO:BP",
                "native": "GO:0000001",
                "name": "test biological process",
                "p_value": 0.01,
                "significant": True,
                "intersection_size": 2,
                "intersections": [
                    [],
                    ["IEA"],
                    ["IMP", "TAS"],
                ],
            },
            {
                "source": "GO:BP",
                "native": "GO:0000002",
                "name": "non-significant process",
                "p_value": 0.08,
                "significant": False,
                "intersection_size": 1,
                "intersections": [
                    ["IEA"],
                    [],
                    [],
                ],
            },
        ],
    }


class TestExtractMembershipRows(unittest.TestCase):

    def test_extracts_members_using_ensg_order(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        self.assertEqual(len(rows), 2)

        self.assertEqual(
            [row["gene_id"] for row in rows],
            ["AT2G00020", "AT3G00030"],
        )

        self.assertEqual(
            [row["evidence_codes"] for row in rows],
            ["IEA", "IMP;TAS"],
        )

    def test_preserves_go_term_and_adjusted_p_value(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        self.assertTrue(
            all(
                row["GO_ID"] == "GO:0000001"
                for row in rows
            )
        )

        self.assertTrue(
            all(
                row["GO_term"] == "test biological process"
                for row in rows
            )
        )

        self.assertTrue(
            all(
                row["adjusted_p_value"] == 0.01
                for row in rows
            )
        )

    def test_excludes_non_significant_terms(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        self.assertNotIn(
            "GO:0000002",
            [row["GO_ID"] for row in rows],
        )

    def test_rejects_intersection_length_mismatch(self):
        response = make_response()
        response["result"][0]["intersections"].pop()

        with self.assertRaisesRegex(
            ValueError,
            "intersection length mismatch",
        ):
            extract_membership_rows(
                response,
                "Conserved 45",
            )

    def test_rejects_intersection_size_mismatch(self):
        response = make_response()
        response["result"][0]["intersection_size"] = 3

        with self.assertRaisesRegex(
            ValueError,
            "intersection_size mismatch",
        ):
            extract_membership_rows(
                response,
                "Conserved 45",
            )

    def test_rejects_unexpected_organism(self):
        response = make_response()
        response["meta"]["query_metadata"]["organism"] = "hsapiens"

        with self.assertRaisesRegex(
            ValueError,
            "unexpected organism",
        ):
            extract_membership_rows(
                response,
                "Conserved 45",
            )

    def test_rejects_unexpected_significance_threshold(self):
        response = make_response()
        response["meta"]["query_metadata"]["user_threshold"] = 0.1

        with self.assertRaisesRegex(
            ValueError,
            "expected significance threshold 0.05",
        ):
            extract_membership_rows(
                response,
                "Conserved 45",
            )


class TestBuildMembershipTable(unittest.TestCase):

    def test_builds_table_from_all_saved_responses(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)

            for _, filename in GENE_SETS:
                json_path = output_dir / filename
                json_path.write_text(
                    json.dumps(make_response()),
                    encoding="utf-8",
                )

            dataframe, summaries = build_membership_table(
                output_dir=output_dir,
            )

        self.assertEqual(
            dataframe.columns.tolist(),
            OUTPUT_COLUMNS,
        )

        self.assertEqual(len(dataframe), 6)
        self.assertEqual(len(summaries), 3)

        self.assertEqual(
            dataframe["gene_set"].nunique(),
            3,
        )

    def test_rejects_duplicate_gene_term_memberships(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)

            response = make_response()
            duplicate_term = response["result"][0].copy()
            response["result"].append(duplicate_term)

            for _, filename in GENE_SETS:
                (output_dir / filename).write_text(
                    json.dumps(response),
                    encoding="utf-8",
                )

            with self.assertRaisesRegex(
                ValueError,
                "Duplicate gene-term memberships",
            ):
                build_membership_table(
                    output_dir=output_dir,
                )


class TestCompareMembershipTables(unittest.TestCase):

    def test_returns_full_match_for_identical_tables(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        expected = pd.DataFrame(
            rows,
            columns=OUTPUT_COLUMNS,
        )

        comparison = compare_membership_tables(
            expected,
            expected.copy(),
        )

        self.assertTrue(comparison["matches"])
        self.assertEqual(len(comparison["missing"]), 0)
        self.assertEqual(len(comparison["extra"]), 0)
        self.assertEqual(
            len(comparison["changed_p_values"]),
            0,
        )

    def test_detects_missing_membership(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        expected = pd.DataFrame(
            rows,
            columns=OUTPUT_COLUMNS,
        )

        existing = expected.iloc[:1].copy()

        comparison = compare_membership_tables(
            expected,
            existing,
        )

        self.assertFalse(comparison["matches"])
        self.assertEqual(len(comparison["missing"]), 1)
        self.assertEqual(len(comparison["extra"]), 0)

    def test_detects_changed_p_value(self):
        rows = extract_membership_rows(
            make_response(),
            "Conserved 45",
        )

        expected = pd.DataFrame(
            rows,
            columns=OUTPUT_COLUMNS,
        )

        existing = expected.copy()
        existing.loc[0, "adjusted_p_value"] = 0.02

        comparison = compare_membership_tables(
            expected,
            existing,
        )

        self.assertFalse(comparison["matches"])
        self.assertEqual(
            len(comparison["changed_p_values"]),
            1,
        )


if __name__ == "__main__":
    unittest.main()
