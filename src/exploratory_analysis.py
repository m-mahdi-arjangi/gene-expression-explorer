from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


# --------------------------------------------------
# Project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

EXPRESSION_PATH = (
    ROOT / "data" / "processed" / "GSE24177_expression.csv"
)

METADATA_PATH = (
    ROOT / "data" / "metadata" / "sample_metadata.csv"
)

FIGURES_DIR = ROOT / "results" / "figures"


# --------------------------------------------------
# Load data
# --------------------------------------------------

def load_data():
    expression = pd.read_csv(
        EXPRESSION_PATH,
        index_col="probe_id"
    )

    metadata = pd.read_csv(
        METADATA_PATH
    ).set_index("sample_id")

    # Keep metadata in the same order as expression columns.
    metadata = metadata.loc[expression.columns]

    return expression, metadata


# --------------------------------------------------
# 1. Expression distribution
# --------------------------------------------------

def plot_expression_distribution(expression):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.boxplot(
        [expression[column].values for column in expression.columns],
        tick_labels=expression.columns,
        showfliers=False
    )

    ax.set_title("Gene Expression Distribution Across Samples")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Expression Value")

    plt.xticks(rotation=45)
    plt.tight_layout()

    output = FIGURES_DIR / "expression_distribution.png"

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# 2. Sample correlation heatmap
# --------------------------------------------------

def plot_sample_correlation(expression, metadata):
    correlation = expression.corr(method="pearson")

    fig, ax = plt.subplots(figsize=(10, 8))

    image = ax.imshow(
        correlation.values,
        vmin=-1,
        vmax=1,
        aspect="auto"
    )

    labels = [
        f"{sample}\n{metadata.loc[sample, 'condition']}"
        for sample in correlation.columns
    ]

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)

    ax.set_title("Correlation Between Samples")

    fig.colorbar(
        image,
        ax=ax,
        label="Pearson correlation"
    )

    plt.tight_layout()

    output = FIGURES_DIR / "sample_correlation.png"

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# 3. Principal Component Analysis
# --------------------------------------------------

def plot_pca(expression, metadata):
    # Select the 1,000 most variable probes.
    variances = expression.var(axis=1)

    top_probes = variances.nlargest(
        min(1000, len(variances))
    ).index

    selected = expression.loc[top_probes]

    # Samples are rows; probes are columns.
    X = selected.T.to_numpy()

    # PCA centers features internally.
    pca = PCA(n_components=2)

    coordinates = pca.fit_transform(X)

    pca_df = pd.DataFrame(
        coordinates,
        index=selected.columns,
        columns=["PC1", "PC2"]
    )

    pca_df = pca_df.join(metadata)

    fig, ax = plt.subplots(figsize=(10, 7))

    for condition, group in pca_df.groupby("condition"):
        ax.scatter(
            group["PC1"],
            group["PC2"],
            label=condition,
            s=90
        )

        for sample, row in group.iterrows():
            ax.annotate(
                sample,
                (row["PC1"], row["PC2"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8
            )

    variance = pca.explained_variance_ratio_ * 100

    ax.set_xlabel(f"PC1 ({variance[0]:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({variance[1]:.2f}% variance)")

    ax.set_title("PCA of Gene Expression Samples")
    ax.legend(title="Condition")

    plt.tight_layout()

    output = FIGURES_DIR / "pca_samples.png"

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output}")

    print("\nExplained variance:")
    print(f"PC1: {variance[0]:.2f}%")
    print(f"PC2: {variance[1]:.2f}%")

    return pca_df


# --------------------------------------------------
# 4. Mean expression per sample
# --------------------------------------------------

def plot_mean_expression(expression, metadata):
    sample_means = expression.mean(axis=0)

    labels = [
        f"{sample}\n{metadata.loc[sample, 'condition']}"
        for sample in sample_means.index
    ]

    fig, ax = plt.subplots(figsize=(12, 6))

    positions = np.arange(len(sample_means))

    ax.bar(
        positions,
        sample_means.values,
    )

    ax.set_xticks(positions)
    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
    )

    ax.set_title("Mean Gene Expression Across Samples")
    ax.set_xlabel("Sample and Condition")
    ax.set_ylabel("Mean Expression Value")

    fig.tight_layout()

    output = FIGURES_DIR / "mean_expression.png"

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output}")


# --------------------------------------------------
# 5. Expression heatmap
# --------------------------------------------------

def plot_expression_heatmap(expression, metadata):
    # Select the 50 most variable probes.
    variances = expression.var(axis=1)

    top_probes = variances.nlargest(
        min(50, len(variances))
    ).index

    selected = expression.loc[top_probes]

    # Standardize expression values within each probe.
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

    sample_labels = [
        f"{sample}\n{metadata.loc[sample, 'condition']}"
        for sample in standardized.columns
    ]

    fig, ax = plt.subplots(
        figsize=(13, 10)
    )

    image = ax.imshow(
        standardized.values,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
    )

    ax.set_xticks(
        np.arange(len(sample_labels))
    )

    ax.set_xticklabels(
        sample_labels,
        rotation=45,
        ha="right",
        fontsize=8,
    )

    ax.set_yticks(
        np.arange(len(top_probes))
    )

    ax.set_yticklabels(
        top_probes,
        fontsize=6,
    )

    ax.set_title(
        "Expression Heatmap of the 50 Most Variable Probes"
    )

    ax.set_xlabel("Samples")
    ax.set_ylabel("Probe IDs")

    fig.colorbar(
        image,
        ax=ax,
        label="Expression Z-score",
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
# Main
# --------------------------------------------------

def main():
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    expression, metadata = load_data()

    print("=== Exploratory Data Analysis ===")
    print(f"Probes: {expression.shape[0]}")
    print(f"Samples: {expression.shape[1]}")

    print("\nMean expression per sample:")
    print(expression.mean().round(3))

    print("\nStandard deviation per sample:")
    print(expression.std().round(3))

    plot_mean_expression(
    expression,
    metadata
)

    plot_expression_heatmap(
        expression,
        metadata
    )

    plot_expression_distribution(
        expression
    )

    plot_sample_correlation(
        expression,
        metadata
    )

    plot_pca(
        expression,
        metadata
    )

    print("\nEDA completed successfully!")


if __name__ == "__main__":
    main()