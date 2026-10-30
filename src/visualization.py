"""Plotting helpers for expression data."""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_mean_expression(means: pd.Series, output_path: str | Path) -> None:
    """Save a bar chart of mean expression per gene."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 5))
    means.sort_values(ascending=False).plot(kind="bar", ax=ax)
    ax.set_title("Mean Expression per Gene")
    ax.set_xlabel("Gene")
    ax.set_ylabel("Mean expression (arbitrary units)")
    ax.tick_params(axis="x", labelrotation=45)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_expression_heatmap(expression: pd.DataFrame, output_path: str | Path) -> None:
    """Save a heatmap of genes by samples without requiring seaborn."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, max(4, len(expression) * 0.45)))
    image = ax.imshow(expression.to_numpy(), aspect="auto", interpolation="nearest")
    ax.set_title("Gene Expression Matrix")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Gene")
    ax.set_xticks(range(len(expression.columns)), labels=expression.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(expression.index)), labels=expression.index)
    fig.colorbar(image, ax=ax, label="Expression (arbitrary units)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
