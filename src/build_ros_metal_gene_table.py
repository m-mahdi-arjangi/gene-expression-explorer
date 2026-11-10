"""Reconstruct the curated ROS/metal gene table from the conserved signature."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    ROOT_DIR
    / "results"
    / "cross_study"
    / "conserved_drought_signature_45_genes.csv"
)
DEFAULT_OUTPUT = (
    ROOT_DIR
    / "results"
    / "cross_study"
    / "GO_enrichment"
    / "core_ROS_metal_genes.csv"
)

# This is a manually curated candidate list recovered from the original workflow.
SELECTED_GENE_IDS = (
    "AT1G08830",
    "AT1G12520",
    "AT1G75040",
    "AT2G28190",
    "AT3G56240",
    "AT4G23670",
    "AT4G25100",
    "AT5G24380",
)

OUTPUT_COLUMNS = (
    "Gene",
    "gene_description",
    "logFC_GSE24177",
    "FDR_GSE24177",
    "logFC_GSE273926",
    "FDR_GSE273926",
    "regulation_GSE24177",
    "regulation_GSE273926",
    "direction_relation",
)
NUMERIC_COLUMNS = (
    "logFC_GSE24177",
    "FDR_GSE24177",
    "logFC_GSE273926",
    "FDR_GSE273926",
)


class WorkflowValidationError(ValueError):
    """Raised when the source table cannot safely produce the expected output."""


def build_ros_metal_gene_table(source: pd.DataFrame) -> pd.DataFrame:
    """Select and validate the eight curated genes from the 45-gene signature."""
    missing_columns = [
        column for column in OUTPUT_COLUMNS if column not in source.columns
    ]
    if missing_columns:
        raise WorkflowValidationError(
            f"Source table is missing required columns: {missing_columns}"
        )

    # Normalize identifiers only in a copy; never mutate the caller's DataFrame.
    normalized_gene_ids = source["Gene"].astype("string").str.strip()
    selected_mask = normalized_gene_ids.isin(SELECTED_GENE_IDS)
    selected = source.loc[selected_mask, list(OUTPUT_COLUMNS)].copy()
    selected["Gene"] = normalized_gene_ids.loc[selected_mask].astype(str)

    duplicated_ids = sorted(
        selected.loc[selected["Gene"].duplicated(keep=False), "Gene"].unique()
    )
    if duplicated_ids:
        raise WorkflowValidationError(
            f"Duplicate selected gene IDs in source table: {duplicated_ids}"
        )

    found_ids = set(selected["Gene"])
    missing_ids = sorted(set(SELECTED_GENE_IDS) - found_ids)
    if missing_ids:
        raise WorkflowValidationError(
            f"Selected genes are missing from source table: {missing_ids}"
        )

    for column in NUMERIC_COLUMNS:
        numeric_values = pd.to_numeric(selected[column], errors="coerce")
        invalid_mask = numeric_values.isna() | ~np.isfinite(numeric_values)
        if invalid_mask.any():
            bad_genes = selected.loc[invalid_mask, "Gene"].tolist()
            raise WorkflowValidationError(
                f"Column '{column}' has missing or non-finite values for: {bad_genes}"
            )
        selected[column] = numeric_values

    for column in ("FDR_GSE24177", "FDR_GSE273926"):
        invalid_fdr = ~selected[column].between(0.0, 1.0, inclusive="both")
        if invalid_fdr.any():
            bad_genes = selected.loc[invalid_fdr, "Gene"].tolist()
            raise WorkflowValidationError(
                f"Column '{column}' contains values outside [0, 1] for: {bad_genes}"
            )

    return (
        selected.loc[:, list(OUTPUT_COLUMNS)]
        .sort_values("Gene")
        .reset_index(drop=True)
    )


def tables_match(expected: pd.DataFrame, observed: pd.DataFrame) -> bool:
    """Return True only when columns, rows, and parsed values match exactly."""
    if list(expected.columns) != list(observed.columns):
        return False
    if len(expected) != len(observed):
        return False

    expected_sorted = (
        expected.sort_values("Gene").reset_index(drop=True)
        if "Gene" in expected.columns
        else expected.reset_index(drop=True)
    )
    observed_sorted = (
        observed.sort_values("Gene").reset_index(drop=True)
        if "Gene" in observed.columns
        else observed.reset_index(drop=True)
    )
    return expected_sorted.equals(observed_sorted)


def _read_csv(path: Path, description: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeError) as exc:
        raise WorkflowValidationError(
            f"Could not read {description} '{path}': {exc}"
        ) from exc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct and validate core_ROS_metal_genes.csv from the "
            "conserved 45-gene cross-study signature."
        )
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=DEFAULT_SOURCE,
        help=f"Source CSV (default: {DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write the reconstructed table; without this flag files are not changed.",
    )
    args = parser.parse_args(argv)

    try:
        source_df = _read_csv(args.source, "Source table")
        expected = build_ros_metal_gene_table(source_df)
    except (FileNotFoundError, WorkflowValidationError) as exc:
        print(f"ERROR: {exc}")
        return 1

    print(f"Source: {args.source}")
    print(f"Selected genes: {len(expected)}")
    print(f"Expected columns: {len(OUTPUT_COLUMNS)}")

    if not args.write:
        if not args.output.is_file():
            print(
                "DRY RUN: reconstructed table is valid, but the output file "
                "does not exist. Use --write to create it."
            )
            return 0

        try:
            observed = _read_csv(args.output, "Existing output")
        except (FileNotFoundError, WorkflowValidationError) as exc:
            print(f"ERROR: {exc}")
            return 1

        if tables_match(expected, observed):
            print("PASS: Existing output matches the reconstructed table exactly.")
            print(f"Output preserved: {args.output}")
            return 0

        print("FAIL: Existing output differs from the reconstructed table.")
        print(f"Expected shape: {expected.shape}")
        print(f"Observed shape: {observed.shape}")
        print(f"Expected columns: {list(expected.columns)}")
        print(f"Observed columns: {list(observed.columns)}")
        print("No files were changed. Review the difference before using --write.")
        return 1

    # An explicit --write is the only path that can create or replace the CSV.
    if args.output.is_file():
        try:
            observed = _read_csv(args.output, "Existing output")
        except WorkflowValidationError as exc:
            print(f"WARNING: {exc}")
        else:
            if tables_match(expected, observed):
                print("PASS: Existing output already matches; no rewrite was needed.")
                print(f"Output preserved: {args.output}")
                return 0

    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        expected.to_csv(args.output, index=False, encoding="utf-8", lineterminator="\n")
    except OSError as exc:
        print(f"ERROR: Could not write output '{args.output}': {exc}")
        return 1

    print(f"WROTE: {args.output}")
    print(f"Rows: {len(expected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
