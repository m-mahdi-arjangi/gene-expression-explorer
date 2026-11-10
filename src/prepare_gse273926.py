
"""Prepare the GSE273926 supplementary expression workbook.

The workbook is preserved as supplied. This script:
- extracts all 14 samples;
- maps source sample columns to GEO GSM accessions;
- records the FASTQ filenames and experimental metadata;
- creates a separate six-sample undigested matrix for the primary comparison.

Important: outputs retain the measurement values as labeled in the source
workbook. They are not treated as raw counts for DESeq2.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook


# --------------------------------------------------
# Project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT_FILE = (
    ROOT
    / "data"
    / "supplementary"
    / "GSE273926"
    / "GSE273926_TenorioBerrio_2024_processed.xlsx"
)

OUTPUT_DIR = ROOT / "data" / "processed" / "GSE273926"
METADATA_DIR = ROOT / "data" / "metadata"

SHEET_NAME = "TenorioBerrio_normalized_counts"


# --------------------------------------------------
# Sample mapping
# Source column, GEO accession, condition,
# preparation, biological replicate
# --------------------------------------------------

SAMPLE_DEFINITIONS = [
    ("Rep1_ww_undigested", "GSM8439682", "control", "undigested", 1),
    ("Rep2_ww_undigested", "GSM8439683", "control", "undigested", 2),
    ("Rep3_ww_undigested", "GSM8439684", "control", "undigested", 3),
    ("Rep1_md_undigested", "GSM8439685", "mild_drought", "undigested", 1),
    ("Rep2_md_undigested", "GSM8439686", "mild_drought", "undigested", 2),
    ("Rep3_md_undigested", "GSM8439687", "mild_drought", "undigested", 3),
    ("Rep1_ww_digested", "GSM8439688", "control", "digested", 1),
    ("Rep2_ww_digested", "GSM8439689", "control", "digested", 2),
    ("Rep1_md_digested", "GSM8439690", "mild_drought", "digested", 1),
    ("Rep2_md_digested", "GSM8439691", "mild_drought", "digested", 2),
    ("Rep1_ww_fixed_digested", "GSM8439692", "control", "fixed_digested", 1),
    ("Rep2_ww_fixed_digested", "GSM8439693", "control", "fixed_digested", 2),
    ("Rep1_md_fixed_digested", "GSM8439694", "mild_drought", "fixed_digested", 1),
    ("Rep2_md_fixed_digested", "GSM8439695", "mild_drought", "fixed_digested", 2),
]


# --------------------------------------------------
# Workbook loading and validation
# --------------------------------------------------

def load_workbook_data():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Supplementary workbook not found: {INPUT_FILE}"
        )

    workbook = load_workbook(
        filename=INPUT_FILE,
        read_only=True,
        data_only=True,
    )

    try:
        if SHEET_NAME not in workbook.sheetnames:
            raise ValueError(
                f"Expected sheet '{SHEET_NAME}' was not found. "
                f"Available sheets: {workbook.sheetnames}"
            )

        sheet = workbook[SHEET_NAME]

        header = next(
            sheet.iter_rows(
                min_row=1,
                max_row=1,
                values_only=True,
            )
        )

        header = list(header)

        expected_source_columns = [
            item[0] for item in SAMPLE_DEFINITIONS
        ]

        if header[0] != "Gene":
            raise ValueError(
                f"Expected first column 'Gene', found: {header[0]}"
            )

        if header[1:] != expected_source_columns:
            raise ValueError(
                "Workbook sample columns do not match the expected "
                "sample mapping.\n"
                f"Found: {header[1:]}\n"
                f"Expected: {expected_source_columns}"
            )

        # The last row contains FASTQ filenames, not gene measurements.
        last_row = list(
            next(
                sheet.iter_rows(
                    min_row=sheet.max_row,
                    max_row=sheet.max_row,
                    values_only=True,
                )
            )
        )

        if last_row[0] != "Geneid":
            raise ValueError(
                "The last workbook row does not contain the expected "
                "'Geneid' FASTQ metadata marker."
            )

        fastq_by_column = {
            source_column: fastq_file
            for source_column, fastq_file in zip(
                header[1:],
                last_row[1:],
            )
        }

        if any(
            not isinstance(filename, str)
            or not filename.strip()
            for filename in fastq_by_column.values()
        ):
            raise ValueError(
                "One or more FASTQ filenames are missing."
            )

        # Exclude the first header row and the final FASTQ metadata row.
        data_rows = list(
            sheet.iter_rows(
                min_row=2,
                max_row=sheet.max_row - 1,
                values_only=True,
            )
        )

    finally:
        workbook.close()

    if not data_rows:
        raise ValueError("No gene expression rows were found.")

    if any(len(row) != len(header) for row in data_rows):
        raise ValueError(
            "At least one data row has an unexpected number of columns."
        )

    expression = pd.DataFrame(
        data_rows,
        columns=header,
    )

    return expression, fastq_by_column


# --------------------------------------------------
# Expression matrix validation
# --------------------------------------------------

def validate_expression(expression):
    gene_ids = expression["Gene"].astype("string").str.strip()

    if gene_ids.isna().any() or gene_ids.eq("").any():
        raise ValueError("Missing or empty gene identifiers found.")

    duplicate_count = int(gene_ids.duplicated().sum())

    if duplicate_count:
        raise ValueError(
            f"Found {duplicate_count} duplicate gene identifiers."
        )

    source_columns = [
        item[0] for item in SAMPLE_DEFINITIONS
    ]

    numeric_data = expression[source_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    if numeric_data.isna().any().any():
        raise ValueError(
            "Missing or non-numeric expression values were found."
        )

    values = numeric_data.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        raise ValueError(
            "Infinite or non-finite expression values were found."
        )

    if (values < 0).any():
        raise ValueError(
            "Negative expression values were found."
        )

    return gene_ids, numeric_data


# --------------------------------------------------
# Create sample metadata
# --------------------------------------------------

def create_metadata(fastq_by_column):
    records = []

    for (
        source_column,
        sample_id,
        condition,
        preparation,
        replicate,
    ) in SAMPLE_DEFINITIONS:

        fastq_filename = fastq_by_column[source_column]

        records.append(
            {
                "sample_id": sample_id,
                "source_column": source_column,
                "condition": condition,
                "preparation": preparation,
                "replicate": replicate,
                "fastq_filename": fastq_filename,
                "selected_for_primary_comparison": (
                    preparation == "undigested"
                ),
                "organism": "Arabidopsis thaliana",
                "geo_series": "GSE273926",
                "measurement_type": (
                    "normalized_counts_as_labeled_in_source_workbook"
                ),
            }
        )

    metadata = pd.DataFrame(records)

    if metadata["sample_id"].duplicated().any():
        raise ValueError("Duplicate GEO sample IDs in metadata.")

    if metadata["fastq_filename"].duplicated().any():
        raise ValueError("Duplicate FASTQ filenames in metadata.")

    return metadata


# --------------------------------------------------
# Save outputs
# --------------------------------------------------

def main():
    print("=== GSE273926 Data Preparation ===")

    print("\nLoading supplementary Excel workbook...")

    expression, fastq_by_column = load_workbook_data()

    print(f"Genes found: {len(expression):,}")
    print(f"Samples found: {len(SAMPLE_DEFINITIONS)}")

    gene_ids, numeric_data = validate_expression(expression)

    metadata = create_metadata(fastq_by_column)

    # Rename source columns to their official GEO accessions.
    source_to_geo = {
        source_column: sample_id
        for source_column, sample_id, *_ in SAMPLE_DEFINITIONS
    }

    all_samples = numeric_data.rename(
        columns=source_to_geo
    )

    all_samples.insert(
        0,
        "gene_id",
        gene_ids.to_numpy(),
    )

    # Select only the undigested samples for the primary comparison.
    primary_sample_ids = metadata.loc[
        metadata["selected_for_primary_comparison"],
        "sample_id",
    ].tolist()

    primary_samples = all_samples[
        ["gene_id"] + primary_sample_ids
    ].copy()

    # Output directories
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    METADATA_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Output files
    all_output = (
        OUTPUT_DIR
        / "GSE273926_normalized_counts_all_samples.csv"
    )

    primary_output = (
        OUTPUT_DIR
        / "GSE273926_normalized_counts_undigested.csv"
    )

    metadata_output = (
        METADATA_DIR
        / "GSE273926_sample_metadata.csv"
    )

    all_samples.to_csv(
        all_output,
        index=False,
    )

    primary_samples.to_csv(
        primary_output,
        index=False,
    )

    metadata.to_csv(
        metadata_output,
        index=False,
    )

    # Summary
    print("\n=== Validation Summary ===")
    print(f"Unique gene IDs: {gene_ids.nunique():,}")
    print("Duplicate gene IDs: 0")
    print("Missing expression values: 0")
    print("Non-finite expression values: 0")
    print("Negative expression values: 0")

    print("\nSamples per condition and preparation:")
    print(
        metadata.groupby(
            ["preparation", "condition"]
        ).size().to_string()
    )

    print("\nPrimary comparison:")
    print("Control samples: 3")
    print("Mild drought samples: 3")
    print("Preparation: undigested")

    print("\nSaved files:")
    print(all_output)
    print(primary_output)
    print(metadata_output)

    print("\nData preparation completed successfully.")


if __name__ == "__main__":
    main()
