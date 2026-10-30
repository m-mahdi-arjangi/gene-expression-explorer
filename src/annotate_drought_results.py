import gzip
from io import StringIO
from pathlib import Path

import pandas as pd

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent

ANNOTATION_FILE = (
    PROJECT_ROOT / "data" / "raw" / "GPL10948_family.soft.gz"
)

RESULTS_FILE = (
    PROJECT_ROOT / "results" / "limma_early_drought_results.csv"
)

OUTPUT_FILE = (
    PROJECT_ROOT / "results" / "limma_early_drought_annotated.csv"
)


def load_platform_annotation(file_path):
    """Load probe annotations from the GEO platform SOFT file."""

    with gzip.open(
        file_path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as file:
        lines = file.readlines()

    start = next(
        (
            i for i, line in enumerate(lines)
            if line.strip() == "!platform_table_begin"
        ),
        None,
    )

    end = next(
        (
            i for i, line in enumerate(lines)
            if line.strip() == "!platform_table_end"
        ),
        None,
    )

    if start is None or end is None:
        raise ValueError(
            "Could not find the platform annotation table."
        )

    annotation = pd.read_csv(
        StringIO("".join(lines[start + 1:end])),
        sep="\t",
        dtype=str,
    )

    required_columns = {
        "ID",
        "ORF",
        "TAIR_ID",
        "Description",
    }

    missing = required_columns - set(annotation.columns)

    if missing:
        raise ValueError(
            f"Missing annotation columns: {sorted(missing)}"
        )

    annotation = annotation.rename(
        columns={
            "ID": "probe_id",
            "ORF": "locus_tag",
            "TAIR_ID": "gene_id",
            "Description": "gene_description",
        }
    )

    annotation = annotation[
        [
            "probe_id",
            "locus_tag",
            "gene_id",
            "gene_description",
        ]
    ]

    annotation = annotation.drop_duplicates(
        subset=["probe_id"]
    )

    return annotation


def main():
    print("Loading GEO platform annotation...")

    annotation = load_platform_annotation(ANNOTATION_FILE)

    print(f"Annotation rows: {len(annotation)}")

    print("\nLoading limma results...")

    results = pd.read_csv(RESULTS_FILE)

    required_columns = {
        "probe_id",
        "logFC",
        "P.Value",
        "adj.P.Val",
    }

    missing = required_columns - set(results.columns)

    if missing:
        raise ValueError(
            f"Missing limma columns: {sorted(missing)}"
        )

    # Exclude technical control probes.
    results = results[
        ~results["probe_id"].str.startswith("AFFX-", na=False)
    ].copy()

    # Match each probe to its platform annotation.
    merged = results.merge(
        annotation,
        on="probe_id",
        how="left",
        validate="many_to_one",
    )

    # Identify the direction of differential expression.
    merged["regulation"] = "not_significant"

    significant = merged["adj.P.Val"] < 0.05

    merged.loc[
        significant & (merged["logFC"] > 0),
        "regulation",
    ] = "upregulated"

    merged.loc[
        significant & (merged["logFC"] < 0),
        "regulation",
    ] = "downregulated"

    # Save the complete annotated results.
    merged.to_csv(OUTPUT_FILE, index=False)

    print("\n=== Annotation Summary ===")
    print(f"Analyzed probes: {len(merged)}")
    print(f"Matched gene IDs: {merged['gene_id'].notna().sum()}")
    print(
        "Significant probes (adjusted p < 0.05):",
        significant.sum(),
    )
    print(
        "Significant upregulated probes:",
        (merged["regulation"] == "upregulated").sum(),
    )
    print(
        "Significant downregulated probes:",
        (merged["regulation"] == "downregulated").sum(),
    )

    print(f"\nSaved annotated results to:\n{OUTPUT_FILE}")

    print("\n=== Top 10 Results ===")

    print(
        merged.sort_values("adj.P.Val")[
            [
                "probe_id",
                "gene_id",
                "gene_description",
                "logFC",
                "P.Value",
                "adj.P.Val",
                "regulation",
            ]
        ].head(10).to_string(index=False)
    )


if __name__ == "__main__":
    main()