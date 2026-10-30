import gzip
from pathlib import Path

import pandas as pd


# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent

RAW_DATA_PATH = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "GSE24177_series_matrix.txt.gz"
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "processed"


def load_geo_expression_data(file_path):
    """
    Load the expression matrix from a GEO Series Matrix file.

    Returns
    -------
    pandas.DataFrame
        Expression values indexed by probe ID.
    """

    if not file_path.exists():
        raise FileNotFoundError(
            f"Dataset not found: {file_path}"
        )

    # Locate the beginning of the expression table.
    with gzip.open(
        file_path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as file:
        for line_number, line in enumerate(file):
            if line.startswith("!series_matrix_table_begin"):
                header_line = line_number + 1
                break
        else:
            raise ValueError(
                "Could not find the expression table."
            )

    # Read the expression matrix.
    expression_data = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip",
        skiprows=header_line,
        nrows=21517,
        index_col=0,
        quotechar='"',
    )

    # Convert expression values to numeric values.
    expression_data = expression_data.apply(
        pd.to_numeric,
        errors="raise",
    )

    expression_data.index.name = "probe_id"

    return expression_data


def main():
    print("=== GEO Expression Data Loader ===\n")

    expression_data = load_geo_expression_data(
        RAW_DATA_PATH
    )

    print(f"Number of probes: {expression_data.shape[0]}")
    print(f"Number of samples: {expression_data.shape[1]}")

    print("\nFirst five rows:")
    print(expression_data.head())

    print("\nMissing values:")
    print(expression_data.isna().sum().sum())

    print("\nDuplicate probe IDs:")
    print(expression_data.index.duplicated().sum())

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = OUTPUT_DIR / "GSE24177_expression.csv"

    expression_data.to_csv(output_path)

    print(f"\nExpression matrix saved to:\n{output_path}")


if __name__ == "__main__":
    main()