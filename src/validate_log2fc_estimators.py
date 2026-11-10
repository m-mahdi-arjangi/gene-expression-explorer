"""Rebuild and validate GSE273926 published and reconstructed log2FC outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT_DIR = Path(__file__).resolve().parents[1]
BASE_DIR = ROOT_DIR / "results" / "GSE273926" / "article_reference"

GROUPS = ("U", "D", "FD")

SAMPLE_COLUMNS = {
    "U": {
        "ww": ["U_WW_1-1", "U_WW_2-1", "U_WW_2-2"],
        "md": ["U_MD_1-1", "U_MD_2-1", "U_MD_2-2"],
    },
    "D": {
        "ww": ["D_WW_1", "D_WW_2"],
        "md": ["D_MD_1", "D_MD_2"],
    },
    "FD": {
        "ww": ["FD_WW_1", "FD_WW_2"],
        "md": ["FD_MD_1", "FD_MD_2"],
    },
}

ESTIMATOR_COLUMNS = [
    "Gene",
    "U_S4_log2FC",
    "U_mean_log2FC",
    "D_S4_log2FC",
    "D_mean_log2FC",
    "FD_S4_log2FC",
    "FD_mean_log2FC",
]

DIRECTION_OUTPUT_COLUMNS = [
    "Gene",
    "U_article_Log2FC",
    "U_S4_log2FC",
    "U_mean_log2FC",
    "D_S4_log2FC",
    "D_mean_log2FC",
    "FD_S4_log2FC",
    "FD_mean_log2FC",
    "U_article_direction",
    "U_S4_direction",
    "U_mean_log_direction",
    "U_direction_relation",
    "D_S4_direction",
    "D_mean_log_direction",
    "D_direction_relation",
    "FD_S4_direction",
    "FD_mean_log_direction",
    "FD_direction_relation",
]

OUTPUT_NAMES = [
    "U_published_vs_estimated_log2FC.csv",
    "D_published_vs_estimated_log2FC.csv",
    "FD_published_vs_estimated_log2FC.csv",
    "log2FC_estimator_validation.csv",
    "U_701_direction_from_tableS4.csv",
    "U_701_direction_reconstruction_summary.csv",
    "D_estimator_direction_differences.csv",
    "FD_estimator_direction_differences.csv",
    "article_table_s5_U_vs_D.csv",
    "article_table_s5_U_vs_FD.csv",
    "article_table_s5_direction_comparison.csv",
    "U_vs_D_opposite_direction_genes.csv",
]


def require_columns(
    dataframe: pd.DataFrame,
    required: list[str],
    label: str,
) -> None:
    missing = [column for column in required if column not in dataframe.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def normalize_gene_ids(values: pd.Series, label: str) -> pd.Series:
    ids = values.astype("string").str.strip()

    if ids.isna().any() or ids.eq("").any():
        raise ValueError(f"{label} contains blank gene IDs.")

    if ids.duplicated().any():
        duplicates = ids[ids.duplicated()].head(5).tolist()
        raise ValueError(f"{label} contains duplicate gene IDs: {duplicates}")

    return ids.astype(str)


def numeric_values(
    dataframe: pd.DataFrame,
    columns: list[str],
    label: str,
) -> pd.DataFrame:
    result = dataframe.copy()

    for column in columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")

        values = result[column].to_numpy(dtype=float)

        if not np.isfinite(values).all():
            raise ValueError(
                f"{label} contains missing or non-finite values in {column}."
            )

    return result


def validate_table_s4(table_s4: pd.DataFrame) -> pd.DataFrame:
    required = ["Gene"]

    for group in GROUPS:
        required.extend(SAMPLE_COLUMNS[group]["ww"])
        required.extend(SAMPLE_COLUMNS[group]["md"])

    require_columns(table_s4, required, "Table S4")

    data = table_s4.copy()
    data["Gene"] = normalize_gene_ids(data["Gene"], "Table S4")

    expression_columns = [column for column in required if column != "Gene"]
    data = numeric_values(data, expression_columns, "Table S4")

    values = data[expression_columns].to_numpy(dtype=float)

    if (values < 0).any():
        raise ValueError("Table S4 contains negative expression values.")

    return data


def validate_table_s5(
    table: pd.DataFrame,
    group: str,
) -> pd.DataFrame:
    required = [
        "Gene",
        "baseMean",
        "Log2-FC",
        "lfcSE",
        "Stat.",
        "P-value",
        "P (adj.)",
    ]

    if group in ("D", "FD"):
        required.append("DE in U samples")

    require_columns(table, required, f"Table S5 {group}")

    data = table.copy()
    data["Gene"] = normalize_gene_ids(data["Gene"], f"Table S5 {group}")

    numeric_columns = [
        "baseMean",
        "Log2-FC",
        "lfcSE",
        "Stat.",
        "P-value",
        "P (adj.)",
    ]
    data = numeric_values(data, numeric_columns, f"Table S5 {group}")

    if (data["baseMean"] < 0).any():
        raise ValueError(f"Table S5 {group} has negative baseMean values.")

    if (data["lfcSE"] <= 0).any():
        raise ValueError(f"Table S5 {group} has non-positive lfcSE values.")

    for column in ("P-value", "P (adj.)"):
        if not data[column].between(0, 1).all():
            raise ValueError(
                f"Table S5 {group} has {column} values outside [0, 1]."
            )

    if group in ("D", "FD"):
        flags = data["DE in U samples"].astype(str).str.strip()
        if not flags.isin(["Yes", "No"]).all():
            raise ValueError(
                f"Table S5 {group} has invalid 'DE in U samples' flags."
            )
        data["DE in U samples"] = flags

    return data


def classify_direction(values) -> np.ndarray:
    """Classify positive, negative and exact-zero values without coercing zero."""
    array = np.asarray(values, dtype=float)

    if not np.isfinite(array).all():
        raise ValueError("Cannot classify missing or non-finite log2FC values.")

    return np.where(
        array > 0,
        "up",
        np.where(array < 0, "down", "zero"),
    )


def classify_relation(reference, estimate) -> np.ndarray:
    """Compare directions; any exact-zero value makes the relation unresolved."""
    reference_direction = classify_direction(reference)
    estimate_direction = classify_direction(estimate)

    unresolved = (
        (reference_direction == "zero")
        | (estimate_direction == "zero")
    )
    same = (
        (reference_direction == estimate_direction)
        & ~unresolved
    )

    return np.where(
        unresolved,
        "unresolved",
        np.where(same, "same_direction", "opposite_direction"),
    )


def calculate_estimators(table_s4: pd.DataFrame) -> pd.DataFrame:
    """Calculate both Table S4 estimators using a pseudocount of one."""
    data = validate_table_s4(table_s4)

    result = pd.DataFrame({"Gene": data["Gene"]})

    for group in GROUPS:
        ww_columns = SAMPLE_COLUMNS[group]["ww"]
        md_columns = SAMPLE_COLUMNS[group]["md"]

        ww = data[ww_columns].astype(float)
        md = data[md_columns].astype(float)

        ww_mean = ww.mean(axis=1)
        md_mean = md.mean(axis=1)

        result[f"{group}_S4_log2FC"] = np.log2(
            (md_mean + 1) / (ww_mean + 1)
        )

        result[f"{group}_mean_log2FC"] = (
            np.log2(md + 1).mean(axis=1)
            - np.log2(ww + 1).mean(axis=1)
        )

    return result[ESTIMATOR_COLUMNS]


def prepare_published_estimates(
    table_s5: pd.DataFrame,
    estimators: pd.DataFrame,
    group: str,
) -> pd.DataFrame:
    data = validate_table_s5(table_s5, group)

    if group in ("D", "FD"):
        data = data.loc[data["DE in U samples"].eq("Yes")].copy()

    selected = data[["Gene", "Log2-FC", "P (adj.)"]]

    merged = selected.merge(
        estimators,
        on="Gene",
        how="inner",
        validate="one_to_one",
        sort=False,
    )

    if len(merged) != len(selected):
        missing = set(selected["Gene"]) - set(merged["Gene"])
        raise ValueError(
            f"{group}: {len(missing)} published genes are missing from Table S4."
        )

    return merged[
        [
            "Gene",
            "Log2-FC",
            "P (adj.)",
            f"{group}_S4_log2FC",
            f"{group}_mean_log2FC",
        ]
    ].copy()


def calculate_metrics(reference, estimate) -> dict:
    reference = np.asarray(reference, dtype=float)
    estimate = np.asarray(estimate, dtype=float)

    valid = np.isfinite(reference) & np.isfinite(estimate)
    reference = reference[valid]
    estimate = estimate[valid]

    if len(reference) < 2:
        raise ValueError("At least two finite gene pairs are required.")

    errors = estimate - reference

    if np.ptp(reference) == 0 or np.ptp(estimate) == 0:
        pearson_r = np.nan
    else:
        pearson_r = float(stats.pearsonr(reference, estimate).statistic)

    if len(np.unique(reference)) < 2 or len(np.unique(estimate)) < 2:
        spearman_rho = np.nan
    else:
        spearman_rho = float(stats.spearmanr(reference, estimate).statistic)

    resolved = (reference != 0) & (estimate != 0)

    if resolved.any():
        direction_agreement = float(
            (np.sign(reference[resolved]) == np.sign(estimate[resolved])).mean()
            * 100
        )
    else:
        direction_agreement = np.nan

    return {
        "genes_compared": int(len(reference)),
        "pearson_r": pearson_r,
        "spearman_rho": spearman_rho,
        "MAE": float(np.mean(np.abs(errors))),
        "RMSE": float(np.sqrt(np.mean(errors ** 2))),
        "mean_bias": float(np.mean(errors)),
        "direction_agreement_pct": direction_agreement,
    }


def build_validation_table(
    published_estimates: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []

    for group in GROUPS:
        data = published_estimates[group]
        reference = data["Log2-FC"].to_numpy(dtype=float)

        for method, column in (
            (1, f"{group}_S4_log2FC"),
            (2, f"{group}_mean_log2FC"),
        ):
            metrics = calculate_metrics(
                reference,
                data[column].to_numpy(dtype=float),
            )
            rows.append({
                "sample_type": group,
                "method": method,
                **metrics,
            })

    return pd.DataFrame(
        rows,
        columns=[
            "sample_type",
            "method",
            "genes_compared",
            "pearson_r",
            "spearman_rho",
            "MAE",
            "RMSE",
            "mean_bias",
            "direction_agreement_pct",
        ],
    )


def build_direction_table(
    published_u: pd.DataFrame,
    estimators: pd.DataFrame,
) -> pd.DataFrame:
    u = validate_table_s5(published_u, "U")

    data = u[["Gene", "Log2-FC"]].rename(
        columns={"Log2-FC": "U_article_Log2FC"}
    )

    data = data.merge(
        estimators,
        on="Gene",
        how="inner",
        validate="one_to_one",
        sort=False,
    )

    if len(data) != len(u):
        raise ValueError("Not all U genes could be matched to Table S4.")

    data["U_article_direction"] = classify_direction(
        data["U_article_Log2FC"]
    )

    for group in GROUPS:
        data[f"{group}_S4_direction"] = classify_direction(
            data[f"{group}_S4_log2FC"]
        )
        data[f"{group}_mean_log_direction"] = classify_direction(
            data[f"{group}_mean_log2FC"]
        )
        data[f"{group}_direction_relation"] = classify_relation(
            data["U_article_Log2FC"],
            data[f"{group}_S4_log2FC"],
        )

    return data[DIRECTION_OUTPUT_COLUMNS].copy()


def build_direction_summary(direction_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    article_direction = direction_table["U_article_direction"].to_numpy()

    for group in GROUPS:
        estimate_direction = direction_table[
            f"{group}_S4_direction"
        ].to_numpy()

        unresolved = (
            (article_direction == "zero")
            | (estimate_direction == "zero")
        )
        same = (
            (article_direction == estimate_direction)
            & ~unresolved
        )
        opposite = (
            (article_direction != estimate_direction)
            & ~unresolved
        )

        rows.append({
            "sample_type": group,
            "genes_compared": int(len(direction_table)),
            "same_direction": int(same.sum()),
            "opposite_direction": int(opposite.sum()),
            "unresolved": int(unresolved.sum()),
        })

    return pd.DataFrame(
        rows,
        columns=[
            "sample_type",
            "genes_compared",
            "same_direction",
            "opposite_direction",
            "unresolved",
        ],
    )


def build_estimator_direction_differences(
    direction_table: pd.DataFrame,
    group: str,
) -> pd.DataFrame:
    s4_column = f"{group}_S4_direction"
    mean_column = f"{group}_mean_log_direction"

    different = direction_table[s4_column] != direction_table[mean_column]

    return direction_table.loc[
        different,
        DIRECTION_OUTPUT_COLUMNS,
    ].copy().reset_index(drop=True)


def build_cross_group_table(
    published_u: pd.DataFrame,
    published_other: pd.DataFrame,
    group: str,
) -> pd.DataFrame:
    u = validate_table_s5(published_u, "U")
    other = validate_table_s5(published_other, group)

    other = other.loc[other["DE in U samples"].eq("Yes")].copy()

    merged = u[["Gene", "Log2-FC"]].merge(
        other[
            [
                "Gene",
                "baseMean",
                "Log2-FC",
                "lfcSE",
                "Stat.",
                "P-value",
                "P (adj.)",
                "DE in U samples",
            ]
        ],
        on="Gene",
        how="inner",
        validate="one_to_one",
        sort=False,
        suffixes=("_U", f"_{group}"),
    )

    other_logfc = f"Log2-FC_{group}"

    merged["direction_relation"] = classify_relation(
        merged["Log2-FC_U"],
        merged[other_logfc],
    )

    columns = [
        "Gene",
        "Log2-FC_U",
        "baseMean",
        other_logfc,
        "lfcSE",
        "Stat.",
        "P-value",
        "P (adj.)",
        "DE in U samples",
        "direction_relation",
    ]

    return merged[columns].copy()


def build_cross_group_direction_summary(
    cross_tables: dict[str, pd.DataFrame],
    table_s5: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows = []

    for group in ("D", "FD"):
        comparison = f"U_vs_{group}"
        data = cross_tables[group]
        total_other = len(table_s5[group])

        same = int(data["direction_relation"].eq("same_direction").sum())
        opposite = int(
            data["direction_relation"].eq("opposite_direction").sum()
        )
        unresolved = int(data["direction_relation"].eq("unresolved").sum())
        shared = len(data)

        rows.append({
            "comparison": comparison,
            "DEGs_in_other_type": int(total_other),
            "shared_genes_with_U": int(shared),
            "same_direction": same,
            "opposite_direction": opposite,
            "zero_or_unresolved": unresolved,
            "shared_pct_of_other_DEGs": (
                shared / total_other * 100 if total_other else np.nan
            ),
            "same_direction_pct_of_shared": (
                same / shared * 100 if shared else np.nan
            ),
        })

    return pd.DataFrame(
        rows,
        columns=[
            "comparison",
            "DEGs_in_other_type",
            "shared_genes_with_U",
            "same_direction",
            "opposite_direction",
            "zero_or_unresolved",
            "shared_pct_of_other_DEGs",
            "same_direction_pct_of_shared",
        ],
    )


def build_opposite_u_d_genes(
    published_u: pd.DataFrame,
    published_d: pd.DataFrame,
) -> pd.DataFrame:
    u = validate_table_s5(published_u, "U")
    d = validate_table_s5(published_d, "D")
    d = d.loc[d["DE in U samples"].eq("Yes")].copy()

    merged = u[
        ["Gene", "Log2-FC", "P (adj.)"]
    ].merge(
        d[["Gene", "Log2-FC", "P (adj.)"]],
        on="Gene",
        how="inner",
        validate="one_to_one",
        suffixes=("_U", "_D"),
        sort=False,
    )

    relation = classify_relation(
        merged["Log2-FC_U"],
        merged["Log2-FC_D"],
    )

    result = merged.loc[
        relation == "opposite_direction",
        [
            "Gene",
            "Log2-FC_U",
            "Log2-FC_D",
            "P (adj.)_U",
            "P (adj.)_D",
        ],
    ].copy()

    return result.reset_index(drop=True)


def build_outputs(
    table_s4: pd.DataFrame,
    table_s5_u: pd.DataFrame,
    table_s5_d: pd.DataFrame,
    table_s5_fd: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    table_s4 = validate_table_s4(table_s4)

    table_s5 = {
        "U": validate_table_s5(table_s5_u, "U"),
        "D": validate_table_s5(table_s5_d, "D"),
        "FD": validate_table_s5(table_s5_fd, "FD"),
    }

    estimators = calculate_estimators(table_s4)

    published_estimates = {
        group: prepare_published_estimates(
            table_s5[group],
            estimators,
            group,
        )
        for group in GROUPS
    }

    validation = build_validation_table(published_estimates)

    direction_table = build_direction_table(
        table_s5["U"],
        estimators,
    )

    direction_summary = build_direction_summary(direction_table)

    cross_tables = {
        group: build_cross_group_table(
            table_s5["U"],
            table_s5[group],
            group,
        )
        for group in ("D", "FD")
    }

    direction_comparison = build_cross_group_direction_summary(
        cross_tables,
        table_s5,
    )

    opposite_u_d = build_opposite_u_d_genes(
        table_s5["U"],
        table_s5["D"],
    )

    outputs = {
        "U_published_vs_estimated_log2FC.csv": published_estimates["U"],
        "D_published_vs_estimated_log2FC.csv": published_estimates["D"],
        "FD_published_vs_estimated_log2FC.csv": published_estimates["FD"],
        "log2FC_estimator_validation.csv": validation,
        "U_701_direction_from_tableS4.csv": direction_table,
        "U_701_direction_reconstruction_summary.csv": direction_summary,
        "D_estimator_direction_differences.csv": (
            build_estimator_direction_differences(direction_table, "D")
        ),
        "FD_estimator_direction_differences.csv": (
            build_estimator_direction_differences(direction_table, "FD")
        ),
        "article_table_s5_U_vs_D.csv": cross_tables["D"],
        "article_table_s5_U_vs_FD.csv": cross_tables["FD"],
        "article_table_s5_direction_comparison.csv": direction_comparison,
        "U_vs_D_opposite_direction_genes.csv": opposite_u_d,
    }

    return outputs


def dataframes_equal(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    rtol: float = 1e-9,
    atol: float = 1e-10,
) -> bool:
    try:
        pd.testing.assert_frame_equal(
            left.reset_index(drop=True),
            right.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            check_like=False,
            rtol=rtol,
            atol=atol,
        )
        return True
    except AssertionError:
        return False



def dataframes_equal_by_gene(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    rtol: float = 1e-9,
    atol: float = 1e-10,
) -> bool:
    """Compare tables by unique Gene IDs while ignoring row order."""
    if list(left.columns) != list(right.columns):
        return False

    if "Gene" not in left.columns:
        return False

    left = left.copy()
    right = right.copy()

    left["Gene"] = left["Gene"].astype(str)
    right["Gene"] = right["Gene"].astype(str)

    if left["Gene"].duplicated().any():
        return False

    if right["Gene"].duplicated().any():
        return False

    if set(left["Gene"]) != set(right["Gene"]):
        return False

    left = left.sort_values("Gene", kind="mergesort").reset_index(drop=True)
    right = right.sort_values("Gene", kind="mergesort").reset_index(drop=True)

    return dataframes_equal(
        left,
        right,
        rtol=rtol,
        atol=atol,
    )


def matches_expected_scientific_correction(
    filename: str,
    existing: pd.DataFrame,
    rebuilt: pd.DataFrame,
) -> bool:
    """Allow only the three documented zero-direction corrections."""
    corrected = existing.copy()

    if filename == "U_701_direction_from_tableS4.csv":
        mask = corrected["Gene"].astype(str).eq("AT2G16740")
        if int(mask.sum()) != 1:
            return False
        if corrected.loc[mask, "FD_direction_relation"].iloc[0] != "opposite_direction":
            return False
        corrected.loc[mask, "FD_direction_relation"] = "unresolved"
        return dataframes_equal(corrected, rebuilt)

    if filename == "U_701_direction_reconstruction_summary.csv":
        mask = corrected["sample_type"].astype(str).eq("FD")
        if int(mask.sum()) != 1:
            return False

        old = corrected.loc[mask].iloc[0]
        if not (
            int(old["same_direction"]) == 661
            and int(old["opposite_direction"]) == 40
            and int(old["unresolved"]) == 0
        ):
            return False

        corrected.loc[mask, "opposite_direction"] = 39
        corrected.loc[mask, "unresolved"] = 1
        return dataframes_equal(corrected, rebuilt)

    if filename == "FD_estimator_direction_differences.csv":
        mask = corrected["Gene"].astype(str).eq("AT2G16740")
        if int(mask.sum()) != 1:
            return False
        if corrected.loc[mask, "FD_direction_relation"].iloc[0] != "opposite_direction":
            return False
        corrected.loc[mask, "FD_direction_relation"] = "unresolved"
        return dataframes_equal(corrected, rebuilt)

    return False


def compare_existing_outputs(
    outputs: dict[str, pd.DataFrame],
    output_dir: Path,
) -> tuple[list[str], list[str]]:
    needs_update = []
    conflicts = []

    print("\n=== Existing output comparison ===")

    for filename in OUTPUT_NAMES:
        rebuilt = outputs[filename]
        path = output_dir / filename

        if not path.exists():
            print(f"MISSING - new output: {filename}")
            needs_update.append(filename)
            continue

        existing = pd.read_csv(path)

        if filename == "U_vs_D_opposite_direction_genes.csv":
            if dataframes_equal_by_gene(existing, rebuilt):
                print(f"FULL MATCH (gene-keyed; row order ignored): {filename}")
                continue
        elif dataframes_equal(existing, rebuilt):
            print(f"FULL MATCH: {filename}")
            continue

        if matches_expected_scientific_correction(
            filename,
            existing,
            rebuilt,
        ):
            print(f"EXPECTED ZERO-DIRECTION CORRECTION: {filename}")
            needs_update.append(filename)
            continue

        print(f"UNEXPECTED DIFFERENCE: {filename}")
        print("  Existing shape:", existing.shape)
        print("  Rebuilt shape:", rebuilt.shape)
        if list(existing.columns) != list(rebuilt.columns):
            print("  Existing columns:", existing.columns.tolist())
            print("  Rebuilt columns:", rebuilt.columns.tolist())
        conflicts.append(filename)

    return needs_update, conflicts


def load_source_tables() -> tuple[pd.DataFrame, ...]:
    paths = {
        "S4": BASE_DIR / "article_table_s4.csv",
        "U": BASE_DIR / "article_table_s5_U.csv",
        "D": BASE_DIR / "article_table_s5_D.csv",
        "FD": BASE_DIR / "article_table_s5_FD.csv",
    }

    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required supplementary files are missing:\n"
            + "\n".join(missing)
        )

    return (
        pd.read_csv(paths["S4"]),
        pd.read_csv(paths["U"]),
        pd.read_csv(paths["D"]),
        pd.read_csv(paths["FD"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild and validate published and Table S4-derived "
            "log2FC outputs for GSE273926."
        )
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help=(
            "Save missing outputs and the explicitly documented "
            "zero-direction corrections after comparison."
        ),
    )
    args = parser.parse_args()

    table_s4, table_s5_u, table_s5_d, table_s5_fd = load_source_tables()

    outputs = build_outputs(
        table_s4,
        table_s5_u,
        table_s5_d,
        table_s5_fd,
    )

    print("=== Rebuilt Log2FC validation summary ===")
    print(
        outputs["log2FC_estimator_validation.csv"].to_string(index=False)
    )

    print("\n=== Rebuilt 701-gene direction summary ===")
    print(
        outputs["U_701_direction_reconstruction_summary.csv"].to_string(
            index=False
        )
    )

    print("\n=== Rebuilt estimator-direction difference counts ===")
    print(
        "D:",
        len(outputs["D_estimator_direction_differences.csv"]),
    )
    print(
        "FD:",
        len(outputs["FD_estimator_direction_differences.csv"]),
    )

    needs_update, conflicts = compare_existing_outputs(
        outputs,
        BASE_DIR,
    )

    if conflicts:
        print("\nNo files were written because unexpected differences exist:")
        for filename in conflicts:
            print(" -", filename)
        raise SystemExit(1)

    if not args.write:
        print("\nDry run completed. No files were modified.")
        if needs_update:
            print("Files to create or correct with --write:")
            for filename in needs_update:
                print(" -", filename)
        else:
            print("All outputs already match the rebuilt results.")
        return

    if not needs_update:
        print("\nAll outputs are already up to date. Nothing to write.")
        return

    for filename in needs_update:
        output_path = BASE_DIR / filename
        outputs[filename].to_csv(output_path, index=False)
        print(f"SAVED: {output_path}")

    print("\nAll requested Log2FC outputs were saved successfully.")


if __name__ == "__main__":
    main()
