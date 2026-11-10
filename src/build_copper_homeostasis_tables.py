"""Reconstruct the exploratory copper-homeostasis comparison tables.

The 15-gene candidate list and sample assignments reproduce the original
manual workflow recovered from the project's PowerShell history.

Default mode validates the existing CSV outputs and never writes files.
Use --write to create missing outputs. Existing mismatched outputs are
protected unless --force is supplied together with --write.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]

LIMMA_INPUT = ROOT_DIR / "results" / "limma_early_drought_annotated.csv"
EXPRESSION_INPUT = (
    ROOT_DIR
    / "data"
    / "processed"
    / "GSE273926"
    / "GSE273926_normalized_counts_undigested.csv"
)
S5_INPUT = (
    ROOT_DIR
    / "results"
    / "GSE273926"
    / "article_reference"
    / "article_table_s5_U.csv"
)
OUTPUT_DIR = ROOT_DIR / "results" / "cross_study" / "GO_enrichment"
REGULON_OUTPUT = OUTPUT_DIR / "copper_homeostasis_regulon_comparison.csv"
SAMPLE_LEVEL_OUTPUT = OUTPUT_DIR / "copper_homeostasis_sample_level.csv"

# Curated candidate panel recovered from the original exploratory workflow.
CANDIDATES: tuple[tuple[str, str], ...] = (
    ("AT5G18830", "SPL7"),
    ("AT5G11260", "HY5"),
    ("AT1G71200", "CITF1"),
    ("AT4G28790", "bHLH23"),
    ("AT5G59030", "COPT1"),
    ("AT3G46900", "COPT2"),
    ("AT2G26975", "COPT6"),
    ("AT5G23980", "FRO4"),
    ("AT5G23990", "FRO5"),
    ("AT1G08830", "CSD1"),
    ("AT2G28190", "CSD2"),
    ("AT1G12520", "CCS1"),
    ("AT4G25100", "FSD1"),
    ("AT3G56240", "CCH"),
    ("AT5G24380", "YSL2"),
)

CONTROL_SAMPLES: tuple[str, ...] = (
    "GSM8439682",
    "GSM8439683",
    "GSM8439684",
)
DROUGHT_SAMPLES: tuple[str, ...] = (
    "GSM8439685",
    "GSM8439686",
    "GSM8439687",
)
ALL_SAMPLES = CONTROL_SAMPLES + DROUGHT_SAMPLES

REGULON_COLUMNS: tuple[str, ...] = (
    "gene_id",
    "symbol",
    "GSE24177_logFC",
    "GSE24177_FDR",
    "GSE24177_regulation",
    "GSE273926_control_mean",
    "GSE273926_drought_mean",
    "GSE273926_descriptive_ratio",
    "GSE273926_published_log2FC",
    "GSE273926_published_FDR",
    "GSE273926_in_published_DEG_list",
)
SAMPLE_LEVEL_COLUMNS: tuple[str, ...] = (
    "gene_id",
    "symbol",
    *ALL_SAMPLES,
    "article_baseMean",
    "article_log2FC",
    "article_FDR",
)

REGULON_NUMERIC_COLUMNS = (
    "GSE24177_logFC",
    "GSE24177_FDR",
    "GSE273926_control_mean",
    "GSE273926_drought_mean",
    "GSE273926_descriptive_ratio",
    "GSE273926_published_log2FC",
    "GSE273926_published_FDR",
)
SAMPLE_LEVEL_NUMERIC_COLUMNS = (*ALL_SAMPLES, "article_baseMean", "article_log2FC", "article_FDR")


def _require_columns(frame: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _candidate_ids() -> tuple[str, ...]:
    return tuple(gene_id for gene_id, _ in CANDIDATES)


def _find_expression_id_column(frame: pd.DataFrame) -> str:
    for candidate in ("Geneid", "gene_id", "Gene", "gene"):
        if candidate in frame.columns:
            return candidate
    if len(frame.columns) == 0:
        raise ValueError("Expression matrix has no columns.")
    return str(frame.columns[0])


def _index_candidate_rows(
    frame: pd.DataFrame,
    identifiers: pd.Series,
    label: str,
    *,
    require_all: bool = True,
) -> dict[str, pd.Series]:
    candidate_ids = set(_candidate_ids())
    subset = frame.loc[identifiers.isin(candidate_ids)].copy()
    subset_ids = identifiers.loc[subset.index]

    duplicated = subset_ids[subset_ids.duplicated(keep=False)].unique().tolist()
    if duplicated:
        raise ValueError(f"{label} contain duplicate candidate gene IDs: {sorted(duplicated)}")

    rows = {str(gene_id): subset.loc[index] for index, gene_id in subset_ids.items()}
    missing = sorted(candidate_ids - set(rows))
    if require_all and missing:
        raise ValueError(f"{label} is missing candidate genes: {missing}")
    return rows


def _as_finite_numeric(values: pd.Series, label: str, *, nonnegative: bool = False) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    array = numeric.to_numpy(dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(f"{label} contains missing or non-numeric values.")
    if nonnegative and (array < 0).any():
        raise ValueError(f"{label} contains negative expression values.")
    return numeric


def build_copper_homeostasis_tables(
    limma: pd.DataFrame,
    expression: pd.DataFrame,
    table_s5: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the regulon-comparison and sample-level tables in memory.

    Missing candidate entries in Limma or published Table S5 are represented
    as missing values, matching the original exploratory table semantics.
    All 15 candidate genes must be uniquely present in the expression matrix.
    """
    _require_columns(limma, ("gene_id", "logFC", "adj.P.Val", "regulation"), "Limma results")
    _require_columns(table_s5, ("Gene", "baseMean", "Log2-FC", "P (adj.)"), "Published Table S5 U")
    _require_columns(expression, ALL_SAMPLES, "GSE273926 normalized expression matrix")

    limma_ids = limma["gene_id"].astype("string").str.strip()
    limma_rows = _index_candidate_rows(
        limma, limma_ids, "Limma results", require_all=False
    )

    expression_id_column = _find_expression_id_column(expression)
    expression_ids = (
        expression[expression_id_column]
        .astype("string")
        .str.extract(r"(AT[1-5]G[0-9]{5})", expand=False)
    )
    expression_rows = _index_candidate_rows(expression, expression_ids, "Expression matrix")

    s5_ids = table_s5["Gene"].astype("string").str.strip()
    s5_subset = table_s5.loc[s5_ids.isin(_candidate_ids())].copy()
    s5_subset_ids = s5_ids.loc[s5_subset.index]
    duplicate_s5 = s5_subset_ids[s5_subset_ids.duplicated(keep=False)].unique().tolist()
    if duplicate_s5:
        raise ValueError(f"Published Table S5 U contains duplicate candidate gene IDs: {sorted(duplicate_s5)}")
    s5_rows = {
        str(gene_id): s5_subset.loc[index]
        for index, gene_id in s5_subset_ids.items()
    }

    regulon_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []

    for gene_id, symbol in CANDIDATES:
        limma_row = limma_rows.get(gene_id)
        expression_row = expression_rows[gene_id]
        article_row = s5_rows.get(gene_id)

        sample_values = _as_finite_numeric(
            expression_row.loc[list(ALL_SAMPLES)],
            f"Expression values for {gene_id}",
            nonnegative=True,
        )
        control_mean = float(sample_values.loc[list(CONTROL_SAMPLES)].mean())
        drought_mean = float(sample_values.loc[list(DROUGHT_SAMPLES)].mean())
        descriptive_ratio = drought_mean / control_mean if control_mean > 0 else np.nan

        if limma_row is None:
            limma_logfc = np.nan
            limma_fdr = np.nan
            limma_regulation = "not_found"
        else:
            numeric_limma = _as_finite_numeric(
                pd.Series([limma_row["logFC"], limma_row["adj.P.Val"]]),
                f"Limma statistics for {gene_id}",
            )
            limma_logfc = float(numeric_limma.iloc[0])
            limma_fdr = float(numeric_limma.iloc[1])
            limma_regulation = str(limma_row["regulation"])
            if limma_regulation.strip() == "" or limma_regulation.lower() == "nan":
                raise ValueError(f"Limma regulation is missing for {gene_id}.")

        if article_row is None:
            article_base_mean = np.nan
            article_log2fc = np.nan
            article_fdr = np.nan
            in_published_deg_list = False
        else:
            numeric_article = _as_finite_numeric(
                pd.Series([article_row["baseMean"], article_row["Log2-FC"], article_row["P (adj.)"]]),
                f"Published Table S5 statistics for {gene_id}",
            )
            article_base_mean = float(numeric_article.iloc[0])
            article_log2fc = float(numeric_article.iloc[1])
            article_fdr = float(numeric_article.iloc[2])
            in_published_deg_list = True

        regulon_rows.append(
            {
                "gene_id": gene_id,
                "symbol": symbol,
                "GSE24177_logFC": limma_logfc,
                "GSE24177_FDR": limma_fdr,
                "GSE24177_regulation": limma_regulation,
                "GSE273926_control_mean": control_mean,
                "GSE273926_drought_mean": drought_mean,
                "GSE273926_descriptive_ratio": descriptive_ratio,
                "GSE273926_published_log2FC": article_log2fc,
                "GSE273926_published_FDR": article_fdr,
                "GSE273926_in_published_DEG_list": in_published_deg_list,
            }
        )

        sample_row: dict[str, object] = {"gene_id": gene_id, "symbol": symbol}
        sample_row.update({sample: float(sample_values[sample]) for sample in ALL_SAMPLES})
        sample_row.update(
            {
                "article_baseMean": article_base_mean,
                "article_log2FC": article_log2fc,
                "article_FDR": article_fdr,
            }
        )
        sample_rows.append(sample_row)

    regulon = pd.DataFrame(regulon_rows, columns=REGULON_COLUMNS)
    sample_level = pd.DataFrame(sample_rows, columns=SAMPLE_LEVEL_COLUMNS)
    return regulon, sample_level


def _normalized_text(series: pd.Series) -> list[str]:
    return ["" if pd.isna(value) else str(value) for value in series.tolist()]


def _normalized_bool(series: pd.Series) -> list[bool | None]:
    normalized: list[bool | None] = []
    for value in series.tolist():
        if pd.isna(value):
            normalized.append(None)
        elif isinstance(value, (bool, np.bool_)):
            normalized.append(bool(value))
        elif str(value).strip().lower() in {"true", "1"}:
            normalized.append(True)
        elif str(value).strip().lower() in {"false", "0"}:
            normalized.append(False)
        else:
            normalized.append(None)
    return normalized


def tables_match(expected: pd.DataFrame, observed: pd.DataFrame, numeric_columns: Iterable[str]) -> bool:
    """Compare tables with strict column/row order and tolerant float parsing."""
    if list(expected.columns) != list(observed.columns) or len(expected) != len(observed):
        return False

    numeric_set = set(numeric_columns)
    for column in expected.columns:
        left = expected[column]
        right = observed[column]
        if column in numeric_set:
            left_num = pd.to_numeric(left, errors="coerce").to_numpy(dtype=float)
            right_num = pd.to_numeric(right, errors="coerce").to_numpy(dtype=float)
            if not np.allclose(left_num, right_num, rtol=1e-12, atol=1e-12, equal_nan=True):
                return False
        elif column == "GSE273926_in_published_DEG_list":
            if _normalized_bool(left) != _normalized_bool(right):
                return False
        elif _normalized_text(left) != _normalized_text(right):
            return False
    return True


def _read_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in (LIMMA_INPUT, EXPRESSION_INPUT, S5_INPUT):
        if not path.exists():
            raise FileNotFoundError(f"Required input file not found: {path}")
    return (
        pd.read_csv(LIMMA_INPUT),
        pd.read_csv(EXPRESSION_INPUT),
        pd.read_csv(S5_INPUT),
    )


def _handle_output(
    path: Path,
    expected: pd.DataFrame,
    numeric_columns: Iterable[str],
    label: str,
    *,
    write: bool,
    force: bool,
) -> bool:
    if path.exists():
        observed = pd.read_csv(path)
        if tables_match(expected, observed, numeric_columns):
            print(f"PASS: Existing {label} matches the reconstructed table exactly.")
            print(f"Output preserved: {path}")
            return True
        if not write:
            print(f"FAIL: Existing {label} differs from the reconstructed table.")
            print("No file was changed. Review the differences before deciding whether to write.")
            return False
        if not force:
            print(f"REFUSED: {label} exists but differs from the reconstructed table.")
            print("Use --write --force only after reviewing and approving the differences.")
            return False
    elif not write:
        print(f"MISSING: {label} does not exist. Use --write to create it.")
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    expected.to_csv(path, index=False, encoding="utf-8-sig")
    print(f"WROTE: {label}: {path}")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Create missing outputs; matching existing outputs are preserved.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting mismatched outputs (requires --write).",
    )
    args = parser.parse_args(argv)

    if args.force and not args.write:
        parser.error("--force requires --write")

    limma, expression, table_s5 = _read_inputs()
    regulon, sample_level = build_copper_homeostasis_tables(limma, expression, table_s5)

    print(f"Curated candidate genes: {len(CANDIDATES)}")
    print(f"Regulon table columns: {len(regulon.columns)}")
    print(f"Sample-level table columns: {len(sample_level.columns)}")

    results = (
        _handle_output(
            REGULON_OUTPUT,
            regulon,
            REGULON_NUMERIC_COLUMNS,
            "copper_homeostasis_regulon_comparison.csv",
            write=args.write,
            force=args.force,
        ),
        _handle_output(
            SAMPLE_LEVEL_OUTPUT,
            sample_level,
            SAMPLE_LEVEL_NUMERIC_COLUMNS,
            "copper_homeostasis_sample_level.csv",
            write=args.write,
            force=args.force,
        ),
    )
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
