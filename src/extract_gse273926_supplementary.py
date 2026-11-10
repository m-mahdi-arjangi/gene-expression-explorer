"""Extract and validate GSE273926 supplementary Tables S4 and S5."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]

SUPPLEMENTARY_DIR = (
    ROOT / "data" / "supplementary" / "GSE273926"
)

METADATA_DIR = ROOT / "data" / "metadata"

EXPRESSION_FILE = (
    ROOT / "data" / "processed" / "GSE273926"
    / "GSE273926_normalized_counts_all_samples.csv"
)

SAMPLE_METADATA_FILE = (
    METADATA_DIR / "GSE273926_sample_metadata.csv"
)

MEMBERSHIP_FILE = (
    METADATA_DIR / "GSE273926_leaf_transcriptome_membership.csv"
)

S4_FILE = SUPPLEMENTARY_DIR / "NPH-246-840-s010.xlsx"
S5_FILE = SUPPLEMENTARY_DIR / "NPH-246-840-s001.xlsx"

OUTPUT_DIR = (
    ROOT / "results" / "GSE273926" / "article_reference"
)

EXPECTED_S4_ROWS = 32833

EXPECTED_S5_ROWS = {
    "U": 701,
    "D": 2990,
    "FD": 286,
}

S4_SHEET = "table_S4_normalizedbulk_transcr"
S4_HEADER_ROW = 11

S4_COLUMNS = [
    "Gene",
    "U_WW_1-1",
    "U_WW_2-1",
    "U_WW_2-2",
    "U_MD_1-1",
    "U_MD_2-1",
    "U_MD_2-2",
    "D_WW_1",
    "D_WW_2",
    "D_MD_1",
    "D_MD_2",
    "FD_WW_1",
    "FD_WW_2",
    "FD_MD_1",
    "FD_MD_2",
    "U-WW transcriptome?",
]

S5_DE_COLUMNS = [
    "Gene",
    "baseMean",
    "Log2-FC",
    "lfcSE",
    "Stat.",
    "P-value",
    "P (adj.)",
]

S5_COLUMNS = {
    "U": S5_DE_COLUMNS,
    "D": S5_DE_COLUMNS + ["DE in U samples"],
    "FD": S5_DE_COLUMNS + ["DE in U samples"],
}

S5_HEADER_ROWS = {
    "U": 10,
    "D": 5,
    "FD": 5,
}


def require_columns(
    dataframe: pd.DataFrame,
    required: set[str],
    source_name: str,
) -> None:
    """Ensure that a table contains all required columns."""

    missing = required - set(dataframe.columns)

    if missing:
        raise ValueError(
            f"{source_name} is missing columns: "
            f"{sorted(missing)}"
        )


def normalize_gene_ids(series: pd.Series) -> pd.Series:
    """Normalize identifiers without restricting valid ID formats."""

    values = series.dropna().astype(str).str.strip()

    invalid = {"", "nan", "none", "<na>"}

    return values[
        ~values.str.lower().isin(invalid)
    ]


def parse_boolean(value, source_name: str) -> bool:
    """Convert supported spreadsheet Boolean representations."""

    if isinstance(value, (bool, np.bool_)):
        return bool(value)

    if value in (0, 1):
        return bool(value)

    if isinstance(value, str):
        normalized = value.strip().lower()

        if normalized in {"true", "yes", "1"}:
            return True

        if normalized in {"false", "no", "0"}:
            return False

    raise ValueError(
        f"Invalid Boolean value in {source_name}: {value!r}"
    )


def read_table_s4() -> tuple[pd.DataFrame, list[str]]:
    """Extract Table S4 and its final FASTQ metadata row."""

    if not S4_FILE.exists():
        raise FileNotFoundError(
            f"Table S4 workbook not found: {S4_FILE}"
        )

    workbook = load_workbook(
        S4_FILE,
        read_only=True,
        data_only=True,
    )

    try:
        if S4_SHEET not in workbook.sheetnames:
            raise ValueError(
                f"Missing Table S4 sheet: {S4_SHEET}"
            )

        sheet = workbook[S4_SHEET]

        header = list(
            next(
                sheet.iter_rows(
                    min_row=S4_HEADER_ROW,
                    max_row=S4_HEADER_ROW,
                    values_only=True,
                )
            )
        )

        # The first header cell is blank in the original workbook.
        if header[0] is not None:
            raise ValueError(
                "Unexpected first Table S4 header cell. "
                f"Observed: {header[0]!r}"
            )

        header[0] = "Gene"

        if header != S4_COLUMNS:
            raise ValueError(
                "Unexpected Table S4 headers.\n"
                f"Observed: {header}"
            )

        last_row = list(
            next(
                sheet.iter_rows(
                    min_row=sheet.max_row,
                    max_row=sheet.max_row,
                    values_only=True,
                )
            )
        )

        # The final row stores FASTQ filenames, not gene measurements.
        if not last_row or last_row[0] != "Geneid":
            raise ValueError(
                "The final Table S4 row does not contain "
                "the expected Geneid marker."
            )

        fastq_filenames = last_row[1:15]

        if len(fastq_filenames) != 14:
            raise ValueError(
                "Table S4 must contain 14 FASTQ filenames "
                "in its final metadata row."
            )

        if any(
            not isinstance(filename, str)
            or not filename.strip()
            for filename in fastq_filenames
        ):
            raise ValueError(
                "One or more Table S4 FASTQ filenames are missing."
            )

        rows = list(
            sheet.iter_rows(
                min_row=S4_HEADER_ROW + 1,
                max_row=sheet.max_row - 1,
                values_only=True,
            )
        )

    finally:
        workbook.close()

    dataframe = pd.DataFrame(
        rows,
        columns=header,
    )

    if len(dataframe) != EXPECTED_S4_ROWS:
        raise ValueError(
            "Unexpected Table S4 row count: "
            f"{len(dataframe)}; expected {EXPECTED_S4_ROWS}."
        )

    ids = normalize_gene_ids(dataframe["Gene"])

    if len(ids) != len(dataframe):
        raise ValueError(
            "Table S4 contains missing or blank gene identifiers."
        )

    if ids.duplicated().any():
        raise ValueError(
            "Duplicate gene identifiers found in Table S4."
        )

    dataframe["Gene"] = ids.to_numpy()

    sample_columns = S4_COLUMNS[1:15]

    numeric = dataframe[sample_columns].apply(
        pd.to_numeric,
        errors="raise",
    )

    values = numeric.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        raise ValueError(
            "Table S4 contains non-finite expression values."
        )

    if (values < 0).any():
        raise ValueError(
            "Table S4 contains negative expression values."
        )

    dataframe[sample_columns] = numeric

    membership = dataframe[
        "U-WW transcriptome?"
    ].map(
        lambda value: parse_boolean(
            value,
            "U-WW transcriptome?",
        )
    )

    dataframe["U-WW transcriptome?"] = membership

    return dataframe, fastq_filenames


def read_table_s5() -> dict[str, pd.DataFrame]:
    """Extract the U, D, and FD sheets from Supplementary Table S5."""

    if not S5_FILE.exists():
        raise FileNotFoundError(
            f"Table S5 workbook not found: {S5_FILE}"
        )

    workbook = load_workbook(
        S5_FILE,
        read_only=True,
        data_only=True,
    )

    extracted = {}

    try:
        for sample_type in ("U", "D", "FD"):

            if sample_type not in workbook.sheetnames:
                raise ValueError(
                    f"Missing Table S5 sheet: {sample_type}"
                )

            sheet = workbook[sample_type]
            header_row = S5_HEADER_ROWS[sample_type]

            header = list(
                next(
                    sheet.iter_rows(
                        min_row=header_row,
                        max_row=header_row,
                        values_only=True,
                    )
                )
            )

            expected = S5_COLUMNS[sample_type]

            if header != expected:
                raise ValueError(
                    f"Unexpected Table S5 {sample_type} headers.\n"
                    f"Observed: {header}\n"
                    f"Expected: {expected}"
                )

            rows = list(
                sheet.iter_rows(
                    min_row=header_row + 1,
                    max_row=sheet.max_row,
                    values_only=True,
                )
            )

            dataframe = pd.DataFrame(
                rows,
                columns=header,
            ).dropna(how="all")

            expected_rows = EXPECTED_S5_ROWS[sample_type]

            if len(dataframe) != expected_rows:
                raise ValueError(
                    f"Unexpected Table S5 {sample_type} row count: "
                    f"{len(dataframe)}; expected {expected_rows}."
                )

            ids = normalize_gene_ids(dataframe["Gene"])

            if len(ids) != len(dataframe):
                raise ValueError(
                    f"Blank gene identifiers in Table S5 {sample_type}."
                )

            if ids.duplicated().any():
                raise ValueError(
                    f"Duplicate gene identifiers in Table S5 {sample_type}."
                )

            dataframe["Gene"] = ids.to_numpy()

            numeric_columns = [
                "baseMean",
                "Log2-FC",
                "lfcSE",
                "Stat.",
                "P-value",
                "P (adj.)",
            ]

            numeric = dataframe[numeric_columns].apply(
                pd.to_numeric,
                errors="raise",
            )

            numeric_values = numeric.to_numpy(dtype=float)

            if not np.isfinite(numeric_values).all():
                raise ValueError(
                    f"Non-finite numeric values in Table S5 {sample_type}."
                )

            if (numeric["baseMean"] < 0).any():
                raise ValueError(
                    f"Negative baseMean values in Table S5 {sample_type}."
                )

            if (numeric["lfcSE"] <= 0).any():
                raise ValueError(
                    f"Non-positive lfcSE values in Table S5 {sample_type}."
                )

            for column in ("P-value", "P (adj.)"):
                if not numeric[column].between(0, 1).all():
                    raise ValueError(
                        f"Invalid probability values in "
                        f"Table S5 {sample_type}: {column}."
                    )

            dataframe[numeric_columns] = numeric

            if "DE in U samples" in dataframe.columns:
                flags = (
                    dataframe["DE in U samples"]
                    .astype(str)
                    .str.strip()
                    .str.lower()
                )

                if not flags.isin({"yes", "no"}).all():
                    raise ValueError(
                        f"Invalid DE in U samples values "
                        f"in Table S5 {sample_type}."
                    )

                dataframe["DE in U samples"] = flags.map(
                    {"yes": "Yes", "no": "No"}
                )

            extracted[sample_type] = dataframe

    finally:
        workbook.close()

    return extracted


def validate_membership(table_s4: pd.DataFrame) -> None:
    """Compare Table S4 membership with the existing metadata CSV."""

    if not MEMBERSHIP_FILE.exists():
        raise FileNotFoundError(
            f"Membership file not found: {MEMBERSHIP_FILE}"
        )

    membership = pd.read_csv(MEMBERSHIP_FILE)

    required = {
        "gene_id",
        "expressed_in_undigested_ww_leaf_transcriptome",
    }

    require_columns(
        membership,
        required,
        "Leaf transcriptome membership",
    )

    ids = normalize_gene_ids(membership["gene_id"])

    if len(ids) != len(membership):
        raise ValueError(
            "Blank identifiers found in membership file."
        )

    if ids.duplicated().any():
        raise ValueError(
            "Duplicate identifiers found in membership file."
        )

    if not ids.reset_index(drop=True).equals(
        table_s4["Gene"].reset_index(drop=True)
    ):
        raise ValueError(
            "Membership gene IDs or their ordering "
            "differ from Table S4."
        )

    existing_flags = membership[
        "expressed_in_undigested_ww_leaf_transcriptome"
    ].map(
        lambda value: parse_boolean(
            value,
            "existing membership CSV",
        )
    )

    extracted_flags = table_s4[
        "U-WW transcriptome?"
    ].reset_index(drop=True)

    if not np.array_equal(
        existing_flags.to_numpy(dtype=bool),
        extracted_flags.to_numpy(dtype=bool),
    ):
        raise ValueError(
            "Table S4 membership flags differ from "
            "the existing metadata CSV."
        )


def build_sample_mapping(
    table_s4: pd.DataFrame,
    fastq_filenames: list[str],
) -> pd.DataFrame:
    """Match each Table S4 sample to one GEO GSM by exact values."""

    if not EXPRESSION_FILE.exists():
        raise FileNotFoundError(
            f"Processed expression matrix not found: {EXPRESSION_FILE}"
        )

    if not SAMPLE_METADATA_FILE.exists():
        raise FileNotFoundError(
            f"Sample metadata not found: {SAMPLE_METADATA_FILE}"
        )

    matrix = pd.read_csv(EXPRESSION_FILE)
    metadata = pd.read_csv(SAMPLE_METADATA_FILE)

    require_columns(
        matrix,
        {"gene_id"},
        "Processed expression matrix",
    )

    metadata_required = {
        "sample_id",
        "condition",
        "preparation",
        "replicate",
        "fastq_filename",
        "experiment_round",
    }

    require_columns(
        metadata,
        metadata_required,
        "Sample metadata",
    )

    matrix_ids = normalize_gene_ids(matrix["gene_id"])
    table_s4_ids = normalize_gene_ids(table_s4["Gene"])

    if matrix_ids.duplicated().any():
        raise ValueError(
            "Duplicate gene identifiers in processed expression matrix."
        )

    if not table_s4_ids.reset_index(drop=True).equals(
        matrix_ids.reset_index(drop=True)
    ):
        raise ValueError(
            "Table S4 gene identifiers or ordering do not "
            "match the processed expression matrix."
        )

    gsm_columns = [
        column for column in matrix.columns
        if column != "gene_id"
    ]

    if len(gsm_columns) != 14:
        raise ValueError(
            f"Expected 14 GSM columns; found {len(gsm_columns)}."
        )

    metadata_ids = (
        metadata["sample_id"]
        .astype(str)
        .str.strip()
    )

    if metadata_ids.duplicated().any():
        raise ValueError(
            "Duplicate sample IDs in sample metadata."
        )

    if set(metadata_ids) != set(gsm_columns):
        raise ValueError(
            "Sample metadata GSM identifiers do not match "
            "the expression matrix columns."
        )

    metadata["sample_id"] = metadata_ids

    for column in gsm_columns:
        matrix[column] = pd.to_numeric(
            matrix[column],
            errors="raise",
        )

        if not np.isfinite(
            matrix[column].to_numpy(dtype=float)
        ).all():
            raise ValueError(
                f"Non-finite expression values in {column}."
            )

    table_sample_columns = S4_COLUMNS[1:15]

    records = []
    used_gsms = set()

    for index, s4_column in enumerate(table_sample_columns):

        s4_values = pd.to_numeric(
            table_s4[s4_column],
            errors="raise",
        ).to_numpy(dtype=float)

        matches = []

        for gsm in gsm_columns:

            matrix_values = matrix[gsm].to_numpy(dtype=float)

            if np.array_equal(s4_values, matrix_values):
                matches.append(gsm)

        if len(matches) != 1:
            raise ValueError(
                f"Table S4 sample {s4_column} has "
                f"{len(matches)} exact GSM matches; expected one."
            )

        gsm = matches[0]

        if gsm in used_gsms:
            raise ValueError(
                f"Multiple Table S4 columns map to {gsm}."
            )

        used_gsms.add(gsm)

        records.append({
            "Table_S4_column": s4_column,
            "GSM": gsm,
            "Table_S4_FASTQ": fastq_filenames[index],
        })

    if used_gsms != set(gsm_columns):
        raise ValueError(
            "The Table S4 to GEO mapping is not one-to-one "
            "across all 14 samples."
        )

    mapping = pd.DataFrame(records)

    mapping = mapping.merge(
        metadata[
            [
                "sample_id",
                "condition",
                "preparation",
                "replicate",
                "fastq_filename",
                "experiment_round",
            ]
        ],
        left_on="GSM",
        right_on="sample_id",
        how="left",
        validate="one_to_one",
    )

    if mapping["sample_id"].isna().any():
        raise ValueError(
            "One or more GSM samples have no metadata."
        )

    if not (
        mapping["Table_S4_FASTQ"].astype(str).to_numpy()
        == mapping["fastq_filename"].astype(str).to_numpy()
    ).all():
        raise ValueError(
            "Table S4 FASTQ filenames do not match "
            "the corresponding sample metadata."
        )

    return mapping[
        [
            "Table_S4_column",
            "GSM",
            "condition",
            "preparation",
            "replicate",
            "experiment_round",
            "fastq_filename",
        ]
    ]


def build_s5_summary(
    tables: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Reconstruct the published Table S5 summary."""

    records = []

    for sample_type in ("U", "D", "FD"):
        dataframe = tables[sample_type]

        flags = dataframe.get(
            "DE in U samples",
            pd.Series(dtype="object"),
        )

        records.append({
            "sample_type": sample_type,
            "reported_DEG_rows": len(dataframe),
            "rows_with_FDR_lt_0.05": int(
                (dataframe["P (adj.)"] < 0.05).sum()
            ),
            "positive_log2FC": int(
                (dataframe["Log2-FC"] > 0).sum()
            ),
            "negative_log2FC": int(
                (dataframe["Log2-FC"] < 0).sum()
            ),
            "min_adjusted_p": dataframe["P (adj.)"].min(),
            "max_adjusted_p": dataframe["P (adj.)"].max(),
            "DE_in_U_yes": (
                int((flags == "Yes").sum())
                if len(flags) else np.nan
            ),
            "DE_in_U_no": (
                int((flags == "No").sum())
                if len(flags) else np.nan
            ),
        })

    return pd.DataFrame(records)


def build_outputs():
    """Extract, validate, and assemble all reproducible outputs."""

    for path in (S4_FILE, S5_FILE):
        if not path.exists():
            raise FileNotFoundError(
                f"Supplementary source file not found: {path}"
            )

    print("Extracting Table S4...")

    table_s4, fastq_filenames = read_table_s4()

    print("Validating leaf transcriptome membership...")

    validate_membership(table_s4)

    print("Extracting Table S5...")

    tables_s5 = read_table_s5()

    print("Validating Table S4 sample-to-GEO mapping...")

    sample_mapping = build_sample_mapping(
        table_s4,
        fastq_filenames,
    )

    s5_summary = build_s5_summary(tables_s5)

    outputs = {
        "article_table_s4.csv": table_s4,
        "article_table_s5_U.csv": tables_s5["U"],
        "article_table_s5_D.csv": tables_s5["D"],
        "article_table_s5_FD.csv": tables_s5["FD"],
        "article_table_s5_summary.csv": s5_summary,
        "table_s4_sample_mapping.csv": sample_mapping,
    }

    metrics = {
        "table_s4_genes": len(table_s4),
        "table_s4_samples": len(S4_COLUMNS[1:15]),
        "table_s5_U_genes": len(tables_s5["U"]),
        "table_s5_D_genes": len(tables_s5["D"]),
        "table_s5_FD_genes": len(tables_s5["FD"]),
        "unique_GEO_samples_mapped": sample_mapping["GSM"].nunique(),
        "membership_rows_validated": len(table_s4),
    }

    return outputs, metrics


def compare_dataframes(
    existing: pd.DataFrame,
    rebuilt: pd.DataFrame,
) -> bool:
    """Compare complete tables while allowing CSV float precision."""

    if list(existing.columns) != list(rebuilt.columns):
        return False

    if existing.shape != rebuilt.shape:
        return False

    try:
        pd.testing.assert_frame_equal(
            existing.reset_index(drop=True),
            rebuilt.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=1e-9,
            atol=1e-12,
        )
    except AssertionError:
        return False

    return True


def compare_existing_outputs(
    outputs: dict[str, pd.DataFrame],
) -> list[str]:
    """Compare existing CSVs without modifying them."""

    conflicts = []

    print("\n=== Existing supplementary output comparison ===")

    for filename, rebuilt in outputs.items():

        path = OUTPUT_DIR / filename

        if not path.exists():
            print(f"MISSING - new output: {filename}")
            continue

        existing = pd.read_csv(path)

        if compare_dataframes(existing, rebuilt):
            print(f"FULL MATCH: {filename}")
        else:
            print(f"DIFFERENT: {filename}")
            conflicts.append(filename)

    return conflicts


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract and validate GSE273926 supplementary Tables S4 and S5."
        )
    )

    parser.add_argument(
        "--write",
        action="store_true",
        help=(
            "Save extracted CSV files after checking existing outputs."
        ),
    )

    args = parser.parse_args()

    outputs, metrics = build_outputs()

    print("\n=== Supplementary extraction summary ===")

    for name, value in metrics.items():
        print(f"{name}: {value}")

    print("\n=== Table S5 summary ===")
    print(
        outputs["article_table_s5_summary.csv"].to_string(
            index=False
        )
    )

    conflicts = compare_existing_outputs(outputs)

    if args.write and conflicts:
        raise SystemExit(
            "\nExisting supplementary outputs differ from "
            "the reconstructed tables. No files were overwritten. "
            f"Review these files: {conflicts}"
        )

    if args.write:
        OUTPUT_DIR.mkdir(
            parents=True,
            exist_ok=True,
        )

        for filename, dataframe in outputs.items():
            dataframe.to_csv(
                OUTPUT_DIR / filename,
                index=False,
            )

        print(
            "\nAll supplementary tables were saved successfully."
        )

    else:
        print(
            "\nDRY RUN: no files were changed. "
            "Run with --write after reviewing the comparison."
        )


if __name__ == "__main__":
    main()
