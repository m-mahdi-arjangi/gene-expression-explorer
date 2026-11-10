
"""Unit tests for GSE273926 data preparation."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from openpyxl import Workbook

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

import prepare_gse273926 as preparation


# ==================================================
# 1. Workbook loading tests
# ==================================================

class LoadWorkbookDataTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        self.source_columns = [
            item[0] for item in preparation.SAMPLE_DEFINITIONS
        ]

    def create_workbook(
        self,
        sheet_name=None,
        headers=None,
        gene_rows=None,
        fastq_names=None,
        marker="Geneid",
    ):
        path = Path(self.temp_dir.name) / "test_workbook.xlsx"

        workbook = Workbook()
        sheet = workbook.active
        sheet.title = sheet_name or preparation.SHEET_NAME

        if headers is None:
            headers = ["Gene"] + self.source_columns

        sheet.append(headers)

        if gene_rows is None:
            values = list(range(1, len(self.source_columns) + 1))
            gene_rows = [
                ["AT1G01010"] + values
            ]

        for row in gene_rows:
            sheet.append(row)

        if fastq_names is None:
            fastq_names = [
                f"{column}.fastq.gz"
                for column in self.source_columns
            ]

        sheet.append([marker] + fastq_names)

        workbook.save(path)
        workbook.close()

        return path

    def test_loads_valid_workbook(self):
        path = self.create_workbook()

        with patch.object(preparation, "INPUT_FILE", path):
            expression, fastq_mapping = (
                preparation.load_workbook_data()
            )

        self.assertEqual(expression.shape[0], 1)
        self.assertEqual(
            expression["Gene"].iloc[0],
            "AT1G01010"
        )
        self.assertEqual(
            len(fastq_mapping),
            len(self.source_columns)
        )

        self.assertEqual(
            fastq_mapping[self.source_columns[0]],
            f"{self.source_columns[0]}.fastq.gz"
        )

    def test_rejects_missing_workbook(self):
        missing_path = (
            Path(self.temp_dir.name) / "missing.xlsx"
        )

        with patch.object(
            preparation,
            "INPUT_FILE",
            missing_path
        ):
            with self.assertRaises(FileNotFoundError):
                preparation.load_workbook_data()

    def test_rejects_missing_expected_sheet(self):
        path = self.create_workbook(
            sheet_name="UnexpectedSheet"
        )

        with patch.object(preparation, "INPUT_FILE", path):
            with self.assertRaisesRegex(
                ValueError,
                "Expected sheet"
            ):
                preparation.load_workbook_data()

    def test_rejects_unexpected_sample_headers(self):
        headers = ["Gene"] + self.source_columns
        headers[1] = "UnexpectedSample"

        path = self.create_workbook(headers=headers)

        with patch.object(preparation, "INPUT_FILE", path):
            with self.assertRaisesRegex(
                ValueError,
                "sample columns do not match"
            ):
                preparation.load_workbook_data()

    def test_rejects_missing_geneid_marker(self):
        path = self.create_workbook(
            marker="UnexpectedMarker"
        )

        with patch.object(preparation, "INPUT_FILE", path):
            with self.assertRaisesRegex(
                ValueError,
                "Geneid"
            ):
                preparation.load_workbook_data()

    def test_rejects_missing_fastq_filename(self):
        fastq_names = [
            f"{column}.fastq.gz"
            for column in self.source_columns
        ]

        fastq_names[0] = None

        path = self.create_workbook(
            fastq_names=fastq_names
        )

        with patch.object(preparation, "INPUT_FILE", path):
            with self.assertRaisesRegex(
                ValueError,
                "FASTQ filenames are missing"
            ):
                preparation.load_workbook_data()

    def test_rejects_workbook_without_gene_rows(self):
        path = self.create_workbook(gene_rows=[])

        with patch.object(preparation, "INPUT_FILE", path):
            with self.assertRaisesRegex(
                ValueError,
                "No gene expression rows"
            ):
                preparation.load_workbook_data()


# ==================================================
# 2. Expression validation tests
# ==================================================

class ValidateExpressionTests(unittest.TestCase):

    def setUp(self):
        self.source_columns = [
            item[0] for item in preparation.SAMPLE_DEFINITIONS
        ]

    def make_valid_expression(self):
        data = {
            "Gene": [
                "AT1G01010",
                "AT1G01020",
            ]
        }

        for index, column in enumerate(self.source_columns):
            data[column] = [
                float(index + 1),
                float(index + 2),
            ]

        return pd.DataFrame(data)

    def test_accepts_valid_expression(self):
        expression = self.make_valid_expression()

        gene_ids, numeric_data = (
            preparation.validate_expression(expression)
        )

        self.assertEqual(len(gene_ids), 2)
        self.assertEqual(
            list(gene_ids),
            ["AT1G01010", "AT1G01020"]
        )
        self.assertEqual(
            numeric_data.shape,
            (2, len(self.source_columns))
        )

    def test_rejects_missing_gene_ids(self):
        expression = self.make_valid_expression()
        expression.loc[0, "Gene"] = pd.NA

        with self.assertRaisesRegex(
            ValueError,
            "Missing or empty gene identifiers"
        ):
            preparation.validate_expression(expression)

    def test_rejects_blank_gene_ids(self):
        expression = self.make_valid_expression()
        expression.loc[0, "Gene"] = "   "

        with self.assertRaisesRegex(
            ValueError,
            "Missing or empty gene identifiers"
        ):
            preparation.validate_expression(expression)

    def test_rejects_duplicate_gene_ids(self):
        expression = self.make_valid_expression()
        expression.loc[1, "Gene"] = "AT1G01010"

        with self.assertRaisesRegex(
            ValueError,
            "duplicate gene identifiers"
        ):
            preparation.validate_expression(expression)

    def test_rejects_missing_expression_values(self):
        expression = self.make_valid_expression()
        expression.loc[0, self.source_columns[0]] = np.nan

        with self.assertRaisesRegex(
            ValueError,
            "Missing or non-numeric expression values"
        ):
            preparation.validate_expression(expression)


    def test_rejects_non_numeric_expression_values(self):
        expression = self.make_valid_expression()
        column = self.source_columns[0]

        # Allow mixed values in this test column.
        expression[column] = expression[column].astype("object")
        expression.loc[0, column] = "invalid"

        with self.assertRaisesRegex(
            ValueError,
            "Missing or non-numeric expression values"
        ):
            preparation.validate_expression(expression)


    def test_rejects_infinite_expression_values(self):
        expression = self.make_valid_expression()
        expression.loc[0, self.source_columns[0]] = np.inf

        with self.assertRaisesRegex(
            ValueError,
            "Infinite or non-finite expression values"
        ):
            preparation.validate_expression(expression)

    def test_rejects_negative_expression_values(self):
        expression = self.make_valid_expression()
        expression.loc[0, self.source_columns[0]] = -1

        with self.assertRaisesRegex(
            ValueError,
            "Negative expression values"
        ):
            preparation.validate_expression(expression)


# ==================================================
# 3. Metadata creation tests
# ==================================================

class CreateMetadataTests(unittest.TestCase):

    def make_fastq_mapping(self):
        return {
            item[0]: f"{item[0]}.fastq.gz"
            for item in preparation.SAMPLE_DEFINITIONS
        }

    def test_creates_valid_metadata(self):
        metadata = preparation.create_metadata(
            self.make_fastq_mapping()
        )

        self.assertEqual(
            len(metadata),
            len(preparation.SAMPLE_DEFINITIONS)
        )

        self.assertFalse(
            metadata["sample_id"].duplicated().any()
        )

        self.assertFalse(
            metadata["fastq_filename"].duplicated().any()
        )

        self.assertEqual(
            metadata["selected_for_primary_comparison"].sum(),
            sum(
                item[3] == "undigested"
                for item in preparation.SAMPLE_DEFINITIONS
            )
        )

        self.assertTrue(
            (
                metadata["organism"]
                == "Arabidopsis thaliana"
            ).all()
        )

    def test_rejects_duplicate_sample_ids(self):
        duplicate_definitions = [
            (
                "Run1",
                "GSM000001",
                "control",
                "undigested",
                1,
            ),
            (
                "Run2",
                "GSM000001",
                "mild_drought",
                "undigested",
                2,
            ),
        ]

        mapping = {
            "Run1": "run1.fastq.gz",
            "Run2": "run2.fastq.gz",
        }

        with patch.object(
            preparation,
            "SAMPLE_DEFINITIONS",
            duplicate_definitions
        ):
            with self.assertRaisesRegex(
                ValueError,
                "Duplicate GEO sample IDs"
            ):
                preparation.create_metadata(mapping)

    def test_rejects_duplicate_fastq_filenames(self):
        mapping = {
            item[0]: "same.fastq.gz"
            for item in preparation.SAMPLE_DEFINITIONS
        }

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate FASTQ filenames"
        ):
            preparation.create_metadata(mapping)


if __name__ == "__main__":
    unittest.main()
