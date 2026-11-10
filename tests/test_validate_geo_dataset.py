
"""Unit tests for GEO expression and metadata validation."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

from validate_geo_dataset import (
    load_sample_metadata,
    validate_dataset,
)


class ValidateGeoDatasetTests(unittest.TestCase):

    def make_expression(self):
        return pd.DataFrame(
            [
                [10.0, 20.0],
                [30.0, 40.0],
            ],
            index=["AT1G01010", "AT1G01020"],
            columns=["GSM1", "GSM2"],
        )

    def make_metadata(self):
        return pd.DataFrame({
            "sample_id": ["GSM1", "GSM2"],
            "condition": ["control", "drought"],
            "timepoint": ["day1", "day1"],
            "replicate": [1, 1],
        })

    def test_accepts_valid_dataset(self):
        expression = self.make_expression()
        metadata = self.make_metadata()

        result = validate_dataset(expression, metadata)

        self.assertIsNone(result)

    def test_rejects_mismatched_sample_ids(self):
        expression = self.make_expression()
        metadata = self.make_metadata()

        metadata.loc[1, "sample_id"] = "GSM3"

        with self.assertRaisesRegex(
            ValueError,
            "do not match"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_duplicate_probe_ids(self):
        expression = self.make_expression()
        metadata = self.make_metadata()

        expression.index = [
            "AT1G01010",
            "AT1G01010",
        ]

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate probe IDs"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_duplicate_expression_sample_ids(self):
        expression = self.make_expression()
        expression.columns = ["GSM1", "GSM1"]

        metadata = pd.DataFrame({
            "sample_id": ["GSM1"],
            "condition": ["control"],
            "timepoint": ["day1"],
            "replicate": [1],
        })

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate expression sample IDs"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_missing_expression_values(self):
        expression = self.make_expression()
        metadata = self.make_metadata()

        expression.iloc[0, 0] = float("nan")

        with self.assertRaisesRegex(
            ValueError,
            "Missing expression values"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_infinite_expression_values(self):
        expression = self.make_expression()
        metadata = self.make_metadata()

        expression.iloc[0, 0] = float("inf")

        with self.assertRaisesRegex(
            ValueError,
            "Non-finite expression values"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_duplicate_metadata_sample_ids(self):
        expression = pd.DataFrame(
            [[10.0]],
            index=["AT1G01010"],
            columns=["GSM1"],
        )

        metadata = pd.DataFrame({
            "sample_id": ["GSM1", "GSM1"],
            "condition": ["control", "control"],
            "timepoint": ["day1", "day1"],
            "replicate": [1, 2],
        })

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate metadata sample IDs"
        ):
            validate_dataset(expression, metadata)

    def test_rejects_missing_metadata_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:

            path = Path(temp_dir) / "metadata.csv"

            pd.DataFrame({
                "sample_id": ["GSM1"],
                "condition": ["control"],
            }).to_csv(path, index=False)

            with self.assertRaisesRegex(
                ValueError,
                "Missing metadata columns"
            ):
                load_sample_metadata(path)


if __name__ == "__main__":
    unittest.main()
