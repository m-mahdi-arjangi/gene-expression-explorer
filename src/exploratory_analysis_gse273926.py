
"""Exploratory analysis of GSE273926 undigested leaf samples.

Uses the supplementary normalized expression values for descriptive
visualization only. It does not perform differential expression testing.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


# --------------------------------------------------
# Project paths and configuration
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

EXPRESSION_FILE = (
    ROOT
    / "data"
    / "processed"
    / "GSE273926"
    / "GSE273926_normalized_counts_undigested.csv"
)

METADATA_FILE = (
    ROOT
    / "data"
    / "metadata"
    / "GSE273926_sample_metadata.csv"
)

OUTPUT_DIR = ROOT / "results" / "GSE273926"
FIGURES_DIR = OUTPUT_DIR / "figures"

N_TOP_PCA_GENES = 1000
N_HEATMAP_GENES = 50


# --------------------------------------------------
# Load and validate data
# --------------------------------------------------

def load_data():
    if not EXPRESSION_FILE.exists():
        raise FileNotFoundError(EXPRESSION_FILE)

    if not METADATA_FILE.exists():
        raise FileNotFoundError(METADATA_FILE)

    expression = pd.read_csv(
        EXPRESSION_FILE,
        index_col="gene_id",
    )

    metadata = pd.read_csv(
        METADATA_FILE,
    )

    if expression.index.isna().any():
        raise ValueError("Missing gene identifiers.")

    if expression.index.duplicated().any():
        raise ValueError("Duplicate gene identifiers.")

    metadata = metadata.loc[
        metadata["selected_for_primary_comparison"]
    ].copy()

    metadata = metadata.set_index("sample_id")

    if metadata.index.duplicated().any():
        raise ValueError("Duplicate sample IDs in metadata.")

    if set(expression.columns) != set(metadata.index):
        raise ValueError(
            "Expression columns do not match selected metadata samples."
        )

    # Match metadata ordering to expression matrix columns.
    metadata = metadata.loc[expression.columns]

    expression = expression.apply(
        pd.to_numeric,
        errors="coerce",
    )

    values = expression.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        raise ValueError(
            "Missing or non-finite expression values."
        )

    if (values < 0).any():
        raise ValueError("Negative expression values found.")

    group_sizes = metadata["condition"].value_counts()

    expected_groups = {
        "control": 3,
        "mild_drought": 3,
    }

    if group_sizes.to_dict() != expected_groups:
        raise ValueError(
            f"Unexpected experimental group sizes: {group_sizes.to_dict()}"
        )

    return expression, metadata


# --------------------------------------------------
# Expression distribution
# --------------------------------------------------

def plot_expression_distribution(log_expression):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.boxplot(
        [
            log_expression[sample].values
            for sample in log_expression.columns
        ],
        tick_labels=log_expression.columns,
        showfliers=False,
    )

    ax.set_title(
        "GSE273926 Expression Distribution Across Samples"
    )
    ax.set_xlabel("GEO Sample")
    ax.set_ylabel("log2(Normalized Expression + 1)")

    plt.xticks(rotation=45)
    fig.tight_layout()

    output = FIGURES_DIR / "expression_distribution.png"

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# Sample correlation
# --------------------------------------------------

def plot_sample_correlation(log_expression):
    correlation = log_expression.corr(method="pearson")

    fig, ax = plt.subplots(figsize=(9, 7))

    image = ax.imshow(
        correlation.values,
        vmin=-1,
        vmax=1,
        aspect="auto",
    )

    labels = correlation.columns.tolist()

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(
        labels,
        rotation=90,
        fontsize=8,
    )

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(
        labels,
        fontsize=8,
    )

    ax.set_title("GSE273926 Sample Correlation")

    fig.colorbar(
        image,
        ax=ax,
        label="Pearson correlation",
    )

    fig.tight_layout()

    output = FIGURES_DIR / "sample_correlation.png"

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Saved: {output}")

    return correlation


# --------------------------------------------------
# Principal Component Analysis
# --------------------------------------------------


def plot_pca(log_expression, metadata):
    required_columns = [
        "condition",
        "replicate",
        "experiment_round",
    ]

    missing = [
        col for col in required_columns
        if col not in metadata.columns
    ]

    if missing:
        raise ValueError(
            f"Missing metadata columns for PCA: {missing}"
        )

    gene_variances = log_expression.var(axis=1)

    selected_genes = gene_variances.nlargest(
        min(N_TOP_PCA_GENES, len(gene_variances))
    ).index

    selected = log_expression.loc[selected_genes]

    # Samples are observations; genes are features.
    X = selected.T.to_numpy()

    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(X)

    pca_df = pd.DataFrame(
        coordinates,
        index=selected.columns,
        columns=["PC1", "PC2"],
    )

    # Include experimental round in PCA metadata.
    pca_df = pca_df.join(
        metadata[
            ["condition", "replicate", "experiment_round"]
        ]
    )

    # Save PCA coordinates for reproducibility.
    pca_output = OUTPUT_DIR / "pca_coordinates.csv"

    pca_df.to_csv(
        pca_output,
        index_label="sample_id",
    )

    print(f"Saved PCA coordinates: {pca_output}")

    print("\nPCA coordinates:")
    print(pca_df.round(3).to_string())

    # Color represents condition.
    condition_colors = {
        "control": "tab:blue",
        "mild_drought": "tab:orange",
    }

    # Marker shape represents experimental round.
    round_markers = {
        "experiment_1": "o",
        "experiment_2": "s",
    }

    fig, ax = plt.subplots(figsize=(9, 7))

    for (condition, experiment_round), group in pca_df.groupby(
        ["condition", "experiment_round"],
        sort=True,
    ):
        ax.scatter(
            group["PC1"],
            group["PC2"],
            color=condition_colors[condition],
            marker=round_markers[experiment_round],
            label=f"{condition} / {experiment_round}",
            s=110,
            edgecolors="black",
            linewidths=0.6,
        )

        for sample, row in group.iterrows():
            ax.annotate(
                sample,
                (row["PC1"], row["PC2"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )

    variance = pca.explained_variance_ratio_ * 100

    ax.set_xlabel(f"PC1 ({variance[0]:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({variance[1]:.2f}% variance)")

    ax.set_title(
        "PCA of GSE273926 Undigested Leaf Samples\n"
        "Color: Condition | Marker: Experiment Round"
    )

    ax.legend(
        title="Condition / Experiment Round",
        fontsize=8,
    )

    fig.tight_layout()

    output = FIGURES_DIR / "pca_samples_condition_and_round.png"

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output}")

    print("\nExplained variance:")
    print(f"PC1: {variance[0]:.2f}%")
    print(f"PC2: {variance[1]:.2f}%")

    return pca_df



# --------------------------------------------------
# Expression heatmap
# --------------------------------------------------

def plot_expression_heatmap(log_expression):
    gene_variances = log_expression.var(axis=1)

    selected_genes = gene_variances.nlargest(
        min(N_HEATMAP_GENES, len(gene_variances))
    ).index

    selected = log_expression.loc[selected_genes]

    row_means = selected.mean(axis=1)

    row_std = selected.std(axis=1).replace(
        0,
        np.nan,
    )

    standardized = (
        selected
        .sub(row_means, axis=0)
        .div(row_std, axis=0)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(10, 11))

    image = ax.imshow(
        standardized.values,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
    )

    ax.set_xticks(range(len(standardized.columns)))
    ax.set_xticklabels(
        standardized.columns,
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticks(range(len(selected_genes)))
    ax.set_yticklabels(
        selected_genes,
        fontsize=6,
    )

    ax.set_title(
        "Top 50 Most Variable Genes - GSE273926"
    )

    ax.set_xlabel("GEO Samples")
    ax.set_ylabel("Gene IDs")

    fig.colorbar(
        image,
        ax=ax,
        label="Gene-wise Z-score",
    )

    fig.tight_layout()

    output = FIGURES_DIR / "expression_heatmap.png"

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# Sample quality summary
# --------------------------------------------------

def save_sample_summary(expression, log_expression, metadata):
    summary = pd.DataFrame(index=expression.columns)

    summary["condition"] = metadata["condition"]
    summary["replicate"] = metadata["replicate"]

    summary["mean_expression"] = expression.mean(axis=0)
    summary["median_expression"] = expression.median(axis=0)
    summary["standard_deviation"] = expression.std(axis=0)

    summary["mean_log2_expression"] = log_expression.mean(axis=0)
    summary["median_log2_expression"] = log_expression.median(axis=0)

    summary["zero_fraction"] = (expression == 0).mean(axis=0)

    output = OUTPUT_DIR / "sample_qc_summary.csv"

    summary.to_csv(output, index_label="sample_id")

    print(f"Saved: {output}")

    print("\nSample QC Summary:")
    print(summary.round(3).to_string())


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=== GSE273926 Exploratory Analysis ===")

    expression, metadata = load_data()

    print(f"Genes: {expression.shape[0]:,}")
    print(f"Samples: {expression.shape[1]}")

    print("\nExperimental groups:")
    print(metadata["condition"].value_counts())

    # Log transformation is used for visualization only.
    log_expression = np.log2(expression + 1)

    transformed_output = (
        OUTPUT_DIR
        / "GSE273926_log2_normalized_expression.csv"
    )

    log_expression.to_csv(
        transformed_output,
        index_label="gene_id",
    )

    print(f"\nSaved transformed matrix: {transformed_output}")

    save_sample_summary(
        expression,
        log_expression,
        metadata,
    )

    plot_expression_distribution(log_expression)

    correlation = plot_sample_correlation(log_expression)

    print("\nSample correlations:")
    print(correlation.round(3).to_string())

    plot_pca(
        log_expression,
        metadata,
    )

    plot_expression_heatmap(log_expression)

    print("\nExploratory analysis completed successfully.")
    print("No differential expression test was performed.")


if __name__ == "__main__":
    main()
