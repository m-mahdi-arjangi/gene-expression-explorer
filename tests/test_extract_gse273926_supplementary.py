"""Unit tests for GSE273926 supplementary extraction."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from unittest.mock import patch

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

import extract_gse273926_supplementary as extraction


class SupplementaryExtractionTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

        root = Path(self.temp_dir.name)

        self.s4_path = root / "table_s4.xlsx"
        self.s5_path = root / "table_s5.xlsx"
        self.expression_path = root / "expression.csv"
        self.metadata_path = root / "metadata.csv"
        self.membership_path = root / "membership.csv"
        self.output_dir = root / "outputs"

        patchers = [
            patch.object(extraction, "S4_FILE", self.s4_path),
            patch.object(extraction, "S5_FILE", self.s5_path),
            patch.object(
                extraction, "EXPRESSION_FILE", self.expression_path
            ),
            patch.object(
                extraction, "SAMPLE_METADATA_FILE", self.metadata_path
            ),
            patch.object(
                extraction, "MEMBERSHIP_FILE", self.membership_path
            ),
            patch.object(
                extraction, "OUTPUT_DIR", self.output_dir
            ),
            patch.object(extraction, "EXPECTED_S4_ROWS", 2),
            patch.dict(
                extraction.EXPECTED_S5_ROWS,
                {"U": 2, "D": 2, "FD": 2},
            ),
        ]

        for patcher in patchers:
            patcher.start()
            self.addCleanup(patcher.stop)

        self.gene_ids = ["GENE_A", "GENE_B"]

        self.fastq_names = [
            f"file_{index:02d}.fastq"
            for index in range(1, 15)
        ]

        self.write_s4_workbook()
        self.write_s5_workbook()
        self.write_expression_and_metadata()
        self.write_membership()

    # --------------------------------------------------
    # Synthetic input builders
    # --------------------------------------------------

    def write_s4_workbook(
        self,
        gene_ids=None,
        header=None,
        marker="Geneid",
        negative_expression=False,
        invalid_membership=False,
        sheet_name=None,
    ):
        """Create a small synthetic Table S4 workbook."""

        if gene_ids is None:
            gene_ids = self.gene_ids

        workbook = Workbook()
        sheet = workbook.active
        sheet.title = sheet_name or extraction.S4_SHEET

        for _ in range(extraction.S4_HEADER_ROW - 1):
            sheet.append([None] * 16)

        if header is None:
            header = [None] + extraction.S4_COLUMNS[1:]

        sheet.append(header)

        for row_index, gene_id in enumerate(gene_ids):
            values = [
                sample_index + row_index * 100
                for sample_index in range(1, 15)
            ]

            if negative_expression and row_index == 0:
                values[0] = -1

            membership = row_index == 0

            if invalid_membership and row_index == 0:
                membership = "unknown"

            sheet.append(
                [gene_id] + values + [membership]
            )

        sheet.append(
            [marker] + self.fastq_names + [None]
        )

        workbook.save(self.s4_path)
        workbook.close()

    def write_s5_workbook(
        self,
        include_sheets=("U", "D", "FD"),
        header_override=None,
        rows_override=None,
    ):
        """Create synthetic U, D, and FD Table S5 sheets."""

        rows_by_type = {
            "U": [
                ["GENE_A", 10.0, 1.0, 0.1, 3.0, 0.001, 0.01],
                ["GENE_B", 20.0, -1.0, 0.2, -2.0, 0.02, 0.03],
            ],
            "D": [
                ["GENE_A", 10.0, 0.5, 0.1, 5.0, 0.001, 0.01, "Yes"],
                ["GENE_B", 20.0, -0.5, 0.2, -2.0, 0.02, 0.03, "No"],
            ],
            "FD": [
                ["GENE_A", 10.0, 2.0, 0.1, 5.0, 0.001, 0.01, "No"],
                ["GENE_B", 20.0, 3.0, 0.2, 4.0, 0.02, 0.03, "No"],
            ],
        }

        if rows_override:
            rows_by_type.update(rows_override)

        workbook = Workbook()
        workbook.remove(workbook.active)

        for sample_type in include_sheets:
            sheet = workbook.create_sheet(sample_type)

            header_row = extraction.S5_HEADER_ROWS[sample_type]

            for _ in range(header_row - 1):
                sheet.append([None] * len(
                    extraction.S5_COLUMNS[sample_type]
                ))

            header = extraction.S5_COLUMNS[sample_type].copy()

            if header_override and sample_type in header_override:
                header = header_override[sample_type]

            sheet.append(header)

            for row in rows_by_type[sample_type]:
                sheet.append(row)

        workbook.save(self.s5_path)
        workbook.close()

    def write_expression_and_metadata(self):
        """Create a matrix whose sample vectors match Table S4 exactly."""

        matrix = {
            "gene_id": self.gene_ids,
        }

        metadata_records = []

        for index in range(14):
            gsm = f"GSM{index + 1:05d}"

            values = [
                index + 1,
                index + 101,
            ]

            matrix[gsm] = values

            if index < 6:
                preparation = "undigested"
            elif index < 10:
                preparation = "digested"
            else:
                preparation = "fixed_digested"

            if index in (0, 1, 2, 6, 7, 10, 11):
                condition = "control"
            else:
                condition = "mild_drought"

            replicate = index % 3 + 1

            metadata_records.append({
                "sample_id": gsm,
                "condition": condition,
                "preparation": preparation,
                "replicate": replicate,
                "fastq_filename": self.fastq_names[index],
                "experiment_round": (
                    "experiment_1"
                    if replicate == 1
                    else "experiment_2"
                ),
            })

        pd.DataFrame(matrix).to_csv(
            self.expression_path,
            index=False,
        )

        pd.DataFrame(metadata_records).to_csv(
            self.metadata_path,
            index=False,
        )

    def write_membership(self, flags=None):
        if flags is None:
            flags = [True, False]

        pd.DataFrame({
            "gene_id": self.gene_ids,
            "expressed_in_undigested_ww_leaf_transcriptome": flags,
        }).to_csv(
            self.membership_path,
            index=False,
        )

    # --------------------------------------------------
    # Table S4 tests
    # --------------------------------------------------

    def test_extracts_valid_table_s4(self):
        dataframe, fastq = extraction.read_table_s4()

        self.assertEqual(dataframe.shape, (2, 16))
        self.assertEqual(dataframe["Gene"].tolist(), self.gene_ids)
        self.assertEqual(fastq, self.fastq_names)

    def test_rejects_missing_table_s4_sheet(self):
        self.write_s4_workbook(sheet_name="WrongSheet")

        with self.assertRaisesRegex(
            ValueError,
            "Missing Table S4 sheet",
        ):
            extraction.read_table_s4()

    def test_rejects_unexpected_table_s4_header(self):
        header = [None] + extraction.S4_COLUMNS[1:]
        header[1] = "Unexpected_Sample"

        self.write_s4_workbook(header=header)

        with self.assertRaisesRegex(
            ValueError,
            "Unexpected Table S4 headers",
        ):
            extraction.read_table_s4()

    def test_rejects_missing_geneid_marker(self):
        self.write_s4_workbook(marker="WrongMarker")

        with self.assertRaisesRegex(
            ValueError,
            "Geneid marker",
        ):
            extraction.read_table_s4()

    def test_rejects_duplicate_s4_gene_ids(self):
        self.write_s4_workbook(
            gene_ids=["GENE_A", "GENE_A"]
        )

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate gene identifiers",
        ):
            extraction.read_table_s4()

    def test_rejects_negative_s4_expression(self):
        self.write_s4_workbook(negative_expression=True)

        with self.assertRaisesRegex(
            ValueError,
            "negative expression values",
        ):
            extraction.read_table_s4()

    def test_rejects_invalid_s4_membership_value(self):
        self.write_s4_workbook(invalid_membership=True)

        with self.assertRaisesRegex(
            ValueError,
            "Invalid Boolean value",
        ):
            extraction.read_table_s4()

    # --------------------------------------------------
    # Table S5 tests
    # --------------------------------------------------

    def test_extracts_all_table_s5_sheets(self):
        tables = extraction.read_table_s5()

        self.assertEqual(set(tables), {"U", "D", "FD"})
        self.assertEqual(len(tables["U"]), 2)
        self.assertEqual(len(tables["D"]), 2)
        self.assertEqual(len(tables["FD"]), 2)

        self.assertEqual(
            tables["U"]["Gene"].tolist(),
            self.gene_ids,
        )

    def test_rejects_missing_table_s5_sheet(self):
        self.write_s5_workbook(
            include_sheets=("U", "D")
        )

        with self.assertRaisesRegex(
            ValueError,
            "Missing Table S5 sheet: FD",
        ):
            extraction.read_table_s5()

    def test_rejects_unexpected_table_s5_header(self):
        bad_header = extraction.S5_COLUMNS["U"].copy()
        bad_header[0] = "WrongGeneColumn"

        self.write_s5_workbook(
            header_override={"U": bad_header}
        )

        with self.assertRaisesRegex(
            ValueError,
            "Unexpected Table S5 U headers",
        ):
            extraction.read_table_s5()

    def test_rejects_duplicate_s5_gene_ids(self):
        rows = {
            "U": [
                ["GENE_A", 10, 1, 0.1, 3, 0.001, 0.01],
                ["GENE_A", 20, -1, 0.2, -2, 0.02, 0.03],
            ]
        }

        self.write_s5_workbook(rows_override=rows)

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate gene identifiers",
        ):
            extraction.read_table_s5()

    def test_rejects_invalid_s5_adjusted_p_value(self):
        rows = {
            "U": [
                ["GENE_A", 10, 1, 0.1, 3, 0.001, 1.5],
                ["GENE_B", 20, -1, 0.2, -2, 0.02, 0.03],
            ]
        }

        self.write_s5_workbook(rows_override=rows)

        with self.assertRaisesRegex(
            ValueError,
            "Invalid probability values",
        ):
            extraction.read_table_s5()

    def test_rejects_invalid_de_in_u_flag(self):
        rows = {
            "D": [
                ["GENE_A", 10, 0.5, 0.1, 5, 0.001, 0.01, "Maybe"],
                ["GENE_B", 20, -0.5, 0.2, -2, 0.02, 0.03, "No"],
            ]
        }

        self.write_s5_workbook(rows_override=rows)

        with self.assertRaisesRegex(
            ValueError,
            "Invalid DE in U samples",
        ):
            extraction.read_table_s5()

    # --------------------------------------------------
    # Membership and sample mapping tests
    # --------------------------------------------------

    def test_accepts_matching_leaf_transcriptome_membership(self):
        table_s4, _ = extraction.read_table_s4()

        extraction.validate_membership(table_s4)

    def test_rejects_different_membership_flags(self):
        self.write_membership(flags=[False, True])

        table_s4, _ = extraction.read_table_s4()

        with self.assertRaisesRegex(
            ValueError,
            "membership flags differ",
        ):
            extraction.validate_membership(table_s4)

    def test_builds_exact_one_to_one_sample_mapping(self):
        table_s4, fastq = extraction.read_table_s4()

        mapping = extraction.build_sample_mapping(
            table_s4,
            fastq,
        )

        self.assertEqual(len(mapping), 14)
        self.assertEqual(mapping["GSM"].nunique(), 14)

        self.assertEqual(
            mapping.iloc[0]["Table_S4_column"],
            "U_WW_1-1",
        )
        self.assertEqual(
            mapping.iloc[0]["GSM"],
            "GSM00001",
        )
        self.assertEqual(
            mapping.iloc[0]["condition"],
            "control",
        )

    def test_rejects_sample_expression_mismatch(self):
        matrix = pd.read_csv(self.expression_path)

        matrix["GSM00001"] = [999, 888]

        matrix.to_csv(
            self.expression_path,
            index=False,
        )

        table_s4, fastq = extraction.read_table_s4()

        with self.assertRaisesRegex(
            ValueError,
            "exact GSM matches",
        ):
            extraction.build_sample_mapping(
                table_s4,
                fastq,
            )

    def test_rejects_fastq_metadata_mismatch(self):
        metadata = pd.read_csv(self.metadata_path)

        metadata.loc[0, "fastq_filename"] = "wrong_file.fastq"

        metadata.to_csv(
            self.metadata_path,
            index=False,
        )

        table_s4, fastq = extraction.read_table_s4()

        with self.assertRaisesRegex(
            ValueError,
            "FASTQ filenames do not match",
        ):
            extraction.build_sample_mapping(
                table_s4,
                fastq,
            )

    def test_rejects_duplicate_metadata_sample_ids(self):
        metadata = pd.read_csv(self.metadata_path)

        metadata.loc[1, "sample_id"] = metadata.loc[0, "sample_id"]

        metadata.to_csv(
            self.metadata_path,
            index=False,
        )

        table_s4, fastq = extraction.read_table_s4()

        with self.assertRaisesRegex(
            ValueError,
            "Duplicate sample IDs",
        ):
            extraction.build_sample_mapping(
                table_s4,
                fastq,
            )

    # --------------------------------------------------
    # Summary and output tests
    # --------------------------------------------------

    def test_reconstructs_table_s5_summary(self):
        tables = extraction.read_table_s5()

        summary = extraction.build_s5_summary(tables)

        self.assertEqual(len(summary), 3)

        u = summary.loc[
            summary["sample_type"] == "U"
        ].iloc[0]

        self.assertEqual(u["reported_DEG_rows"], 2)
        self.assertEqual(u["positive_log2FC"], 1)
        self.assertEqual(u["negative_log2FC"], 1)

        d = summary.loc[
            summary["sample_type"] == "D"
        ].iloc[0]

        self.assertEqual(d["DE_in_U_yes"], 1)
        self.assertEqual(d["DE_in_U_no"], 1)

    def test_build_outputs_returns_expected_files(self):
        outputs, metrics = extraction.build_outputs()

        expected_files = {
            "article_table_s4.csv",
            "article_table_s5_U.csv",
            "article_table_s5_D.csv",
            "article_table_s5_FD.csv",
            "article_table_s5_summary.csv",
            "table_s4_sample_mapping.csv",
        }

        self.assertEqual(set(outputs), expected_files)

        self.assertEqual(metrics["table_s4_genes"], 2)
        self.assertEqual(metrics["table_s4_samples"], 14)
        self.assertEqual(metrics["unique_GEO_samples_mapped"], 14)

    def test_dataframe_comparison_accepts_equal_tables(self):
        original = pd.DataFrame({
            "Gene": ["GENE_A", "GENE_B"],
            "value": [1.0, 2.0],
        })

        rebuilt = original.copy()

        self.assertTrue(
            extraction.compare_dataframes(
                original,
                rebuilt,
            )
        )

    def test_dataframe_comparison_detects_changed_values(self):
        original = pd.DataFrame({
            "Gene": ["GENE_A", "GENE_B"],
            "value": [1.0, 2.0],
        })

        rebuilt = original.copy()
        rebuilt.loc[0, "value"] = 999

        self.assertFalse(
            extraction.compare_dataframes(
                original,
                rebuilt,
            )
        )


if __name__ == "__main__":
    unittest.main()
