"""Command-line entry point for the Gene Expression Explorer."""

from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR / "src"))

import pandas as pd
from data_loader import load_expression_data
from visualization import plot_mean_expression, plot_expression_heatmap


def main() -> None:
    data_path = ROOT_DIR / "data" / "sample_gene_expression.csv"
    output_dir = ROOT_DIR / "results" / "figures"

    try:
        expression = load_expression_data(data_path)
    except (FileNotFoundError, ValueError, pd.errors.ParserError) as error:
        print(f"Could not load expression data: {error}", file=sys.stderr)
        raise SystemExit(1) from error

    print("\n=== Gene Expression Explorer ===")
    print(f"Genes: {expression.shape[0]}")
    print(f"Samples: {expression.shape[1]}")
    print("\nExpression matrix:")
    print(expression)

    means = expression.mean(axis=1).sort_values(ascending=False)
    standard_deviations = expression.std(axis=1, ddof=1).fillna(0)

    print("\nMean expression per gene (descending):")
    print(means.round(3).to_string())
    print("\nStandard deviation per gene:")
    print(standard_deviations.round(3).to_string())

    threshold = 20.0
    selected = means[means > threshold]
    print(f"\nGenes with mean expression > {threshold:g}:")
    print(selected.to_string() if not selected.empty else "No genes passed the threshold.")

    minimum = expression.min().min()
    maximum = expression.max().max()
    normalized = (
        expression * 0.0
        if maximum == minimum
        else (expression - minimum) / (maximum - minimum)
    )

    print("\nMin–max normalized expression (first 5 genes):")
    print(normalized.head().round(3).to_string())

    output_dir.mkdir(parents=True, exist_ok=True)
    plot_mean_expression(means, output_dir / "mean_expression.png")
    plot_expression_heatmap(expression, output_dir / "expression_heatmap.png")
    print(f"\nCharts saved to: {output_dir.relative_to(ROOT_DIR)}")


if __name__ == "__main__":
    main()
