"""Reproducible cross-study drought gene overlap analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "results" / "cross_study"

EXPRESSION_FILE = (
    ROOT / "data" / "processed" / "GSE273926"
    / "GSE273926_normalized_counts_all_samples.csv"
)
LIMMA_FILE = ROOT / "results" / "limma_early_drought_annotated.csv"
TABLE_S5_FILE = (
    ROOT / "results" / "GSE273926" / "article_reference"
    / "article_table_s5_U.csv"
)

FDR_THRESHOLD = 0.05


def normalize_ids(series: pd.Series) -> pd.Series:
    """Return non-empty, trimmed gene identifiers."""
    values = series.dropna().astype(str).str.strip()

    invalid = {"", "nan", "none", "<na>"}
    return values[~values.str.lower().isin(invalid)]


def require_columns(
    dataframe: pd.DataFrame,
    required: set[str],
    source_name: str,
) -> None:
    missing = required - set(dataframe.columns)

    if missing:
        raise ValueError(
            f"{source_name} is missing columns: {sorted(missing)}"
        )


def load_inputs():
    """Load and validate the three input datasets."""

    for path in (EXPRESSION_FILE, LIMMA_FILE, TABLE_S5_FILE):
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")

    expression = pd.read_csv(
        EXPRESSION_FILE,
        usecols=["gene_id"],
    )

    limma = pd.read_csv(LIMMA_FILE)
    table_s5 = pd.read_csv(TABLE_S5_FILE)

    require_columns(
        limma,
        {
            "gene_id",
            "probe_id",
            "logFC",
            "adj.P.Val",
            "AveExpr",
            "t",
            "P.Value",
            "B",
            "locus_tag",
            "gene_description",
        },
        "GSE24177 Limma results",
    )

    require_columns(
        table_s5,
        {"Gene", "Log2-FC", "P (adj.)"},
        "GSE273926 Table S5 U",
    )

    expression_ids = normalize_ids(expression["gene_id"])

    if expression_ids.duplicated().any():
        raise ValueError(
            "Duplicate gene IDs found in GSE273926 expression matrix."
        )

    limma["gene_id"] = (
        limma["gene_id"].fillna("").astype(str).str.strip()
    )
    limma["logFC"] = pd.to_numeric(
        limma["logFC"], errors="coerce"
    )
    limma["adj.P.Val"] = pd.to_numeric(
        limma["adj.P.Val"], errors="coerce"
    )

    table_s5["Gene"] = (
        table_s5["Gene"].fillna("").astype(str).str.strip()
    )
    table_s5["Log2-FC"] = pd.to_numeric(
        table_s5["Log2-FC"], errors="coerce"
    )
    table_s5["P (adj.)"] = pd.to_numeric(
        table_s5["P (adj.)"], errors="coerce"
    )

    return expression_ids, limma, table_s5


def calculate_results():
    """Rebuild cross-study overlap and statistical validation."""

    expression_ids, limma, table_s5 = load_inputs()

    # Build the common gene universe from all measured expression
    # identifiers and all annotated Limma gene identifiers.
    limma_all_ids = set(normalize_ids(limma["gene_id"]))
    expression_id_set = set(expression_ids)

    universe = expression_id_set & limma_all_ids

    if not universe:
        raise ValueError("The common gene universe is empty.")

    # Select statistically significant GSE24177 genes.
    sig_limma = limma.loc[
        (limma["gene_id"] != "")
        & (limma["adj.P.Val"] < FDR_THRESHOLD)
    ].copy()

    sig_limma = sig_limma.dropna(
        subset=["logFC", "adj.P.Val"]
    )

    if not np.isfinite(sig_limma["logFC"]).all():
        raise ValueError("Non-finite Limma logFC values detected.")

    if sig_limma["gene_id"].duplicated().any():
        raise ValueError(
            "Duplicate significant gene IDs found in Limma results. "
            "Resolve probe-to-gene mapping before continuing."
        )

    # Select significant genes from the published Table S5 U.
    sig_s5 = table_s5.loc[
        (table_s5["Gene"] != "")
        & (table_s5["P (adj.)"] < FDR_THRESHOLD)
    ].copy()

    sig_s5 = sig_s5.dropna(
        subset=["Log2-FC", "P (adj.)"]
    )

    if not np.isfinite(sig_s5["Log2-FC"]).all():
        raise ValueError("Non-finite Table S5 log2FC values detected.")

    if sig_s5["Gene"].duplicated().any():
        raise ValueError(
            "Duplicate significant gene IDs found in Table S5 U."
        )

    # Restrict both significant sets to the common universe.
    limma_ids = set(sig_limma["gene_id"]) & universe
    s5_ids = set(sig_s5["Gene"]) & universe

    if len(limma_ids) != len(sig_limma):
        raise ValueError(
            "Some significant Limma genes are outside the common universe."
        )

    # Prepare GSE24177 results for one-to-one merging.
    limma_merge = sig_limma[
        [
            "probe_id",
            "gene_id",
            "logFC",
            "AveExpr",
            "t",
            "P.Value",
            "adj.P.Val",
            "B",
            "locus_tag",
            "gene_description",
        ]
    ].copy()

    limma_merge = limma_merge.rename(
        columns={
            "gene_id": "Gene",
            "logFC": "logFC_GSE24177",
            "adj.P.Val": "FDR_GSE24177",
        }
    )

    # Prepare published GSE273926 results.
    s5_merge = sig_s5[
        ["Gene", "Log2-FC", "P (adj.)"]
    ].rename(
        columns={
            "Log2-FC": "logFC_GSE273926",
            "P (adj.)": "FDR_GSE273926",
        }
    )

    shared = limma_merge.merge(
        s5_merge,
        on="Gene",
        how="inner",
        validate="one_to_one",
    )

    shared["regulation_GSE24177"] = np.where(
        shared["logFC_GSE24177"] > 0,
        "upregulated",
        np.where(
            shared["logFC_GSE24177"] < 0,
            "downregulated",
            "unresolved",
        ),
    )

    shared["regulation_GSE273926"] = np.where(
        shared["logFC_GSE273926"] > 0,
        "upregulated",
        np.where(
            shared["logFC_GSE273926"] < 0,
            "downregulated",
            "unresolved",
        ),
    )

    same_direction = (
        (shared["regulation_GSE24177"]
         == shared["regulation_GSE273926"])
        & (shared["regulation_GSE24177"] != "unresolved")
    )

    opposite_direction = (
        (shared["regulation_GSE24177"] != "unresolved")
        & (shared["regulation_GSE273926"] != "unresolved")
        & (
            shared["regulation_GSE24177"]
            != shared["regulation_GSE273926"]
        )
    )

    shared["direction_relation"] = np.select(
        [same_direction, opposite_direction],
        ["same_direction", "opposite_direction"],
        default="unresolved",
    )

    # Preserve the existing output schema.
    shared = shared[
        [
            "probe_id",
            "logFC_GSE24177",
            "AveExpr",
            "t",
            "P.Value",
            "FDR_GSE24177",
            "B",
            "locus_tag",
            "Gene",
            "gene_description",
            "logFC_GSE273926",
            "FDR_GSE273926",
            "direction_relation",
            "regulation_GSE24177",
            "regulation_GSE273926",
        ]
    ].copy()

    shared["regulation"] = np.where(
        shared["logFC_GSE24177"] > 0,
        "upregulated",
        np.where(
            shared["logFC_GSE24177"] < 0,
            "downregulated",
            "unresolved",
        ),
    )

    # Restore the historical column order.
    shared = shared[
        [
            "probe_id",
            "logFC_GSE24177",
            "AveExpr",
            "t",
            "P.Value",
            "FDR_GSE24177",
            "B",
            "locus_tag",
            "Gene",
            "gene_description",
            "regulation",
            "logFC_GSE273926",
            "FDR_GSE273926",
            "direction_relation",
            "regulation_GSE24177",
            "regulation_GSE273926",
        ]
    ]

    shared = shared.sort_values(
        ["FDR_GSE24177", "Gene"]
    ).reset_index(drop=True)

    conserved = shared.loc[
        shared["direction_relation"] == "same_direction"
    ].copy()

    upregulated = shared.loc[
        (shared["regulation_GSE24177"] == "upregulated")
        & (shared["regulation_GSE273926"] == "upregulated")
    ].copy()

    downregulated = shared.loc[
        (shared["regulation_GSE24177"] == "downregulated")
        & (shared["regulation_GSE273926"] == "downregulated")
    ].copy()

    discordant = shared.loc[
        shared["direction_relation"] == "opposite_direction"
    ].copy()

    unresolved_count = int(
        (shared["direction_relation"] == "unresolved").sum()
    )

    same_count = len(conserved)
    opposite_count = len(discordant)
    overlap_count = len(shared)

    same_pct = (
        100 * same_count / overlap_count
        if overlap_count else 0.0
    )

    # Historical overlap summary.
    overlap_summary = pd.DataFrame(
        [
            {
                "significant_GSE24177_with_gene_id": len(sig_limma),
                "significant_GSE273926_U": len(sig_s5),
                "shared_significant_genes": overlap_count,
                "same_direction": same_count,
                "opposite_direction": opposite_count,
                "unresolved": unresolved_count,
                "same_direction_pct_of_overlap": same_pct,
            }
        ]
    )

    # Statistical overlap test within the common universe.
    universe_size = len(universe)
    limma_count = len(limma_ids)
    s5_count = len(s5_ids)

    expected_overlap = (
        limma_count * s5_count / universe_size
    )

    fold_enrichment = (
        overlap_count / expected_overlap
        if expected_overlap else np.nan
    )

    hypergeometric_p = hypergeom.sf(
        overlap_count - 1,
        universe_size,
        limma_count,
        s5_count,
    )

    overlap_table = [
        [
            overlap_count,
            limma_count - overlap_count,
        ],
        [
            s5_count - overlap_count,
            universe_size - limma_count - s5_count
            + overlap_count,
        ],
    ]

    overlap_fisher = fisher_exact(
        overlap_table,
        alternative="two-sided",
    )

    # Directional concordance among shared significant genes.
    up_up = int(
        (
            (shared["regulation_GSE24177"] == "upregulated")
            & (shared["regulation_GSE273926"] == "upregulated")
        ).sum()
    )

    up_down = int(
        (
            (shared["regulation_GSE24177"] == "upregulated")
            & (shared["regulation_GSE273926"] == "downregulated")
        ).sum()
    )

    down_up = int(
        (
            (shared["regulation_GSE24177"] == "downregulated")
            & (shared["regulation_GSE273926"] == "upregulated")
        ).sum()
    )

    down_down = int(
        (
            (shared["regulation_GSE24177"] == "downregulated")
            & (shared["regulation_GSE273926"] == "downregulated")
        ).sum()
    )

    direction_table = [
        [up_up, up_down],
        [down_up, down_down],
    ]

    direction_fisher = fisher_exact(
        direction_table,
        alternative="two-sided",
    )

    statistical_validation = pd.DataFrame(
        [
            {
                "common_tested_universe": universe_size,
                "significant_GSE24177": limma_count,
                "significant_GSE273926": s5_count,
                "observed_overlap": overlap_count,
                "expected_overlap": expected_overlap,
                "fold_enrichment": fold_enrichment,
                "hypergeometric_p": hypergeometric_p,
                "overlap_fisher_odds_ratio": overlap_fisher.statistic,
                "overlap_fisher_p": overlap_fisher.pvalue,
                "up_up": up_up,
                "up_down": up_down,
                "down_up": down_up,
                "down_down": down_down,
                "direction_fisher_odds_ratio": direction_fisher.statistic,
                "direction_fisher_p": direction_fisher.pvalue,
            }
        ]
    )

    universe_df = pd.DataFrame(
        {"gene_id": sorted(universe)}
    )

    outputs = {
        "common_tested_gene_universe.csv": (
            universe_df, "gene_id"
        ),
        "cross_study_overlap_summary.csv": (
            overlap_summary, None
        ),
        "cross_study_statistical_validation.csv": (
            statistical_validation, None
        ),
        "significant_drought_genes_shared.csv": (
            shared, "Gene"
        ),
        "conserved_drought_signature_45_genes.csv": (
            conserved, "Gene"
        ),
        "shared_upregulated_drought_genes_29.csv": (
            upregulated, "Gene"
        ),
        "shared_downregulated_drought_genes_16.csv": (
            downregulated, "Gene"
        ),
        "discordant_drought_genes_7.csv": (
            discordant, "Gene"
        ),
    }

    metrics = {
        "universe": universe_size,
        "limma_significant": len(sig_limma),
        "table_s5_significant_total": len(sig_s5),
        "table_s5_inside_universe": s5_count,
        "table_s5_outside_universe": len(sig_s5) - s5_count,
        "shared": overlap_count,
        "same_direction": same_count,
        "opposite_direction": opposite_count,
        "unresolved": unresolved_count,
    }

    return outputs, metrics


def compare_existing_outputs(outputs):
    """Compare complete rebuilt tables against saved CSVs."""

    mismatches = []

    print("\n=== Full comparison with existing saved outputs ===")

    for filename, (rebuilt, key_column) in outputs.items():
        path = OUTPUT_DIR / filename

        if not path.exists():
            print(f"MISSING: {filename}")
            mismatches.append(filename)
            continue

        old = pd.read_csv(path)

        if list(old.columns) != list(rebuilt.columns):
            print(f"DIFFERENT COLUMNS: {filename}")
            mismatches.append(filename)
            continue

        new = rebuilt.copy()

        # Ignore row ordering for tables with a unique gene key.
        if key_column is not None:
            if (
                key_column not in old.columns
                or old[key_column].duplicated().any()
                or new[key_column].duplicated().any()
            ):
                print(f"INVALID KEY: {filename}")
                mismatches.append(filename)
                continue

            old = old.sort_values(
                key_column, kind="stable"
            ).reset_index(drop=True)

            new = new.sort_values(
                key_column, kind="stable"
            ).reset_index(drop=True)

        if old.shape != new.shape:
            print(
                f"DIFFERENT SHAPE: {filename} | "
                f"old={old.shape}, rebuilt={new.shape}"
            )
            mismatches.append(filename)
            continue

        different_columns = []

        for column in new.columns:
            old_values = old[column]
            new_values = new[column]

            is_numeric = (
                pd.api.types.is_numeric_dtype(old_values)
                or pd.api.types.is_numeric_dtype(new_values)
            )

            if is_numeric:
                a = pd.to_numeric(
                    old_values, errors="coerce"
                ).to_numpy(dtype=float)

                b = pd.to_numeric(
                    new_values, errors="coerce"
                ).to_numpy(dtype=float)

                equal = np.allclose(
                    a,
                    b,
                    rtol=1e-6,
                    atol=1e-12,
                    equal_nan=True,
                )
            else:
                a = (
                    old_values.astype("string")
                    .fillna("<NA>")
                    .astype(str)
                    .to_numpy()
                )

                b = (
                    new_values.astype("string")
                    .fillna("<NA>")
                    .astype(str)
                    .to_numpy()
                )

                equal = np.array_equal(a, b)

            if not equal:
                different_columns.append(column)

        if different_columns:
            print(
                f"DIFFERENT VALUES: {filename} | "
                f"columns={different_columns}"
            )
            mismatches.append(filename)
        else:
            print(f"FULL MATCH: {filename}")

    return mismatches

def main():
    parser = argparse.ArgumentParser(
        description="Rebuild and validate cross-study drought overlap."
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write rebuilt CSV files after comparing existing outputs.",
    )
    args = parser.parse_args()

    outputs, metrics = calculate_results()

    print("\n=== Rebuilt cross-study analysis ===")

    for name, value in metrics.items():
        print(f"{name}: {value}")

    print("\n=== Statistical results ===")
    print(
        outputs["cross_study_statistical_validation.csv"][0]
        .to_string(index=False)
    )

    mismatches = compare_existing_outputs(outputs)

    if args.write and mismatches:
        raise SystemExit(
            "\nExisting outputs differ from rebuilt results. "
            "No files were overwritten. Review the differences first."
        )

    if args.write:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

        for filename, (dataframe, _) in outputs.items():
            dataframe.to_csv(
                OUTPUT_DIR / filename,
                index=False,
            )

        print("\nAll cross-study outputs were saved successfully.")
    else:
        print(
            "\nDRY RUN: no files were changed. "
            "Run with --write only after reviewing the comparison."
        )


if __name__ == "__main__":
    main()
