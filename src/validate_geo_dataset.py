from pathlib import Path

import pandas as pd


# --------------------------------------------------
# Project paths
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

EXPRESSION_PATH = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "GSE24177_expression.csv"
)

METADATA_PATH = (
    PROJECT_ROOT
    / "data"
    / "metadata"
    / "sample_metadata.csv"
)


# --------------------------------------------------
# Load data
# --------------------------------------------------

def load_expression_data(path):
    """Load the gene expression matrix."""

    if not path.exists():
        raise FileNotFoundError(
            f"Expression file not found: {path}"
        )

    expression = pd.read_csv(
        path,
        index_col="probe_id"
    )

    return expression


def load_sample_metadata(path):
    """Load sample annotations."""

    if not path.exists():
        raise FileNotFoundError(
            f"Sample metadata file not found: {path}"
        )

    metadata = pd.read_csv(path)

    required_columns = {
        "sample_id",
        "condition",
        "timepoint",
        "replicate",
    }

    missing_columns = required_columns - set(
        metadata.columns
    )

    if missing_columns:
        raise ValueError(
            f"Missing metadata columns: {missing_columns}"
        )

    return metadata


# --------------------------------------------------
# Validate dataset
# --------------------------------------------------

def validate_dataset(expression, metadata):
    """Validate sample identifiers and expression values."""

    print("\n=== GEO Dataset Validation ===")

    print(f"Number of probes: {expression.shape[0]}")
    print(f"Number of samples: {expression.shape[1]}")

    print(
        "\nMissing expression values:",
        expression.isna().sum().sum()
    )

    print(
        "Duplicate probe IDs:",
        expression.index.duplicated().sum()
    )

    print(
        "Duplicate sample IDs in metadata:",
        metadata["sample_id"].duplicated().sum()
    )

    expression_samples = set(expression.columns)

    metadata_samples = set(metadata["sample_id"])

    missing_metadata = expression_samples - metadata_samples

    extra_metadata = metadata_samples - expression_samples

    print("\nSamples without metadata:", missing_metadata)
    print("Metadata samples absent from expression matrix:", extra_metadata)

    if missing_metadata or extra_metadata:
        raise ValueError(
            "Expression columns and metadata sample IDs do not match."
        )

    if not expression.index.is_unique:
        raise ValueError("Duplicate probe IDs detected.")

    if not expression.columns.is_unique:
        raise ValueError("Duplicate expression sample IDs detected.")

    if expression.isna().any().any():
        raise ValueError("Missing expression values detected.")

    if not expression.map(
        lambda value: pd.notna(value)
        and abs(value) != float("inf")
    ).all().all():
        raise ValueError("Non-finite expression values detected.")

    if metadata["sample_id"].duplicated().any():
        raise ValueError("Duplicate metadata sample IDs detected.")

    print("\nDataset validation passed!")


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    expression = load_expression_data(
        EXPRESSION_PATH
    )

    metadata = load_sample_metadata(
        METADATA_PATH
    )

    validate_dataset(
        expression,
        metadata
    )

    # Arrange metadata in the same order as the matrix.
    metadata = metadata.set_index("sample_id")

    metadata = metadata.loc[expression.columns]

    print("\nSample metadata:")
    print(metadata)

    print("\nSamples per condition and timepoint:")

    print(
        metadata.groupby(
            ["condition", "timepoint"]
        ).size()
    )

    print("\nValidation complete.")


if __name__ == "__main__":
    main()