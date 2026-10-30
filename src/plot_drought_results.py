from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# Project paths
ROOT = Path(__file__).resolve().parent.parent

INPUT_FILE = ROOT / "results" / "limma_early_drought_annotated.csv"
OUTPUT_DIR = ROOT / "results" / "figures"
OUTPUT_FILE = OUTPUT_DIR / "top_drought_genes.png"

SIGNIFICANCE_THRESHOLD = 0.05
TOP_N = 10


def main():
    print("=== Differential Expression Visualization ===")

    # Check input file
    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Annotated results file not found: {INPUT_FILE}"
        )

    # Load annotated limma results
    df = pd.read_csv(INPUT_FILE)

    # Validate required columns
    required_columns = {
        "probe_id",
        "gene_description",
        "logFC",
        "adj.P.Val",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    # Convert statistical columns to numeric
    df["logFC"] = pd.to_numeric(
        df["logFC"],
        errors="coerce",
    )

    df["adj.P.Val"] = pd.to_numeric(
        df["adj.P.Val"],
        errors="coerce",
    )

    # Remove rows with missing statistical values
    df = df.dropna(
        subset=["logFC", "adj.P.Val"]
    ).copy()

    # Exclude technical control probes
    df = df[
        ~df["probe_id"].astype(str).str.startswith("AFFX-")
    ].copy()

    # Keep statistically significant probes
    significant = df[
        df["adj.P.Val"] < SIGNIFICANCE_THRESHOLD
    ].copy()

    if significant.empty:
        raise ValueError(
            "No statistically significant probes found."
        )

    # Select the most increased and decreased probes
    increased = significant.nlargest(
        TOP_N,
        "logFC",
    )

    decreased = significant.nsmallest(
        TOP_N,
        "logFC",
    )

    selected = pd.concat(
        [decreased, increased]
    ).copy()

    # Create readable gene labels
    if "gene_id" in selected.columns:
        gene_ids = selected["gene_id"].fillna(
            selected["probe_id"]
        ).astype(str)
    else:
        gene_ids = selected["probe_id"].astype(str)

    descriptions = (
        selected["gene_description"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.slice(0, 45)
    )

    selected["label"] = [
        f"{gene} — {description} [{probe}]"
        if description
        else f"{gene} [{probe}]"
        for gene, description, probe in zip(
            gene_ids,
            descriptions,
            selected["probe_id"].astype(str),
        )
    ]

    # Sort by log fold change
    selected = selected.sort_values("logFC")

    # Create output directory
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # Generate figure
    height = max(8, len(selected) * 0.48)

    fig, ax = plt.subplots(
        figsize=(14, height)
    )

    ax.barh(
        selected["label"],
        selected["logFC"],
    )

    ax.axvline(
        0,
        linewidth=1,
        color="black",
    )

    ax.set_xlabel(
        "Log2 Fold Change (Drought vs Control)"
    )

    ax.set_ylabel("Gene / Probe")

    ax.set_title(
        "Top 10 Upregulated and Downregulated Probes\n"
        "Moderate Drought vs Control — Day 1"
    )

    fig.tight_layout()

    fig.savefig(
        OUTPUT_FILE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    # Print summary
    print(f"\nAnalyzed probes: {len(df)}")
    print(f"Significant probes: {len(significant)}")
    print(f"Selected probes: {len(selected)}")

    print(f"\nPlot saved to:\n{OUTPUT_FILE}")

    print("\nTop 5 upregulated probes:")

    print(
        increased[
            [
                "probe_id",
                "gene_description",
                "logFC",
                "adj.P.Val",
            ]
        ].head(5).to_string(index=False)
    )

    print("\nTop 5 downregulated probes:")

    print(
        decreased[
            [
                "probe_id",
                "gene_description",
                "logFC",
                "adj.P.Val",
            ]
        ].head(5).to_string(index=False)
    )

    print("\nVisualization completed successfully!")


if __name__ == "__main__":
    main()