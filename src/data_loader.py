"""Input loading and validation for gene-expression tables."""

from pathlib import Path
import numpy as np
import pandas as pd


def load_expression_data(path: str | Path) -> pd.DataFrame:
    """Load a CSV with gene_id plus numeric sample columns.

    Returns a DataFrame indexed by gene ID. Missing or non-numeric expression
    values are rejected rather than silently imputed.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input CSV does not exist: {path}")

    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError("The expression CSV is empty.")
    if "gene_id" not in frame.columns:
        raise ValueError("The CSV must contain a 'gene_id' column.")
    if frame["gene_id"].isna().any() or (frame["gene_id"].astype(str).str.strip() == "").any():
        raise ValueError("Gene IDs must not be missing or blank.")
    if frame["gene_id"].duplicated().any():
        raise ValueError("Gene IDs must be unique.")

    sample_columns = [column for column in frame.columns if column != "gene_id"]
    if not sample_columns:
        raise ValueError("The CSV must contain at least one sample column.")

    values = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    if values.isna().any().any():
        raise ValueError("Expression values must be numeric and must not be missing.")
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError("Expression values must be finite numbers.")

    values.index = frame["gene_id"].astype(str)
    values.index.name = "gene_id"
    return values
