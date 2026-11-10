
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# ==================================================
# 1. Project paths
# ==================================================

ROOT_DIR = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    ROOT_DIR
    / "results"
    / "cross_study"
    / "GO_enrichment"
    / "copper_homeostasis_regulon_comparison.csv"
)

OUTPUT_DIR = (
    ROOT_DIR
    / "results"
    / "cross_study"
    / "GO_enrichment"
    / "figures"
)


# ==================================================
# 2. Genes included in Figure 1
# ==================================================

GENES = [
    {
        "gene_id": "AT1G08830",
        "symbol": "CSD1",
        "module": "Cu/Zn-SOD",
        "shared": True,
    },
    {
        "gene_id": "AT2G28190",
        "symbol": "CSD2",
        "module": "Cu/Zn-SOD",
        "shared": True,
    },
    {
        "gene_id": "AT1G12520",
        "symbol": "CCS1",
        "module": "Cu/Zn-SOD",
        "shared": True,
    },
    {
        "gene_id": "AT4G25100",
        "symbol": "FSD1",
        "module": "Fe-SOD",
        "shared": True,
    },
    {
        "gene_id": "AT3G56240",
        "symbol": "CCH",
        "module": "Metal handling",
        "shared": True,
    },
    {
        "gene_id": "AT5G24380",
        "symbol": "YSL2",
        "module": "Metal handling",
        "shared": True,
    },
    {
        "gene_id": "AT1G71200",
        "symbol": "CITF1",
        "module": "Cu uptake",
        "shared": False,
    },
    {
        "gene_id": "AT3G46900",
        "symbol": "COPT2",
        "module": "Cu uptake",
        "shared": False,
    },
]


# ==================================================
# 3. Load and validate results
# ==================================================

def load_data():

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}"
        )

    df = pd.read_csv(INPUT_FILE)

    required_columns = [
        "gene_id",
        "symbol",
        "GSE24177_logFC",
        "GSE24177_FDR",
        "GSE273926_published_log2FC",
        "GSE273926_published_FDR",
    ]

    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    df["gene_id"] = df["gene_id"].astype(str).str.strip()

    selected_ids = {
        gene["gene_id"]
        for gene in GENES
    }

    selected = df[
        df["gene_id"].isin(selected_ids)
    ].copy()

    if selected["gene_id"].duplicated().any():
        raise ValueError(
            "Duplicated gene IDs found in input data."
        )

    found_ids = set(selected["gene_id"])
    missing_genes = selected_ids - found_ids

    if missing_genes:
        raise ValueError(
            f"Missing genes: {sorted(missing_genes)}"
        )

    numeric_columns = [
        "GSE24177_logFC",
        "GSE24177_FDR",
        "GSE273926_published_log2FC",
        "GSE273926_published_FDR",
    ]

    for column in numeric_columns:
        selected[column] = pd.to_numeric(
            selected[column],
            errors="coerce"
        )

    if selected[numeric_columns].isna().any().any():
        raise ValueError(
            "Missing or non-numeric values found."
        )

    selected = selected.set_index("gene_id")

    return selected


# ==================================================
# 4. Plot cross-study log2FC
# ==================================================

def plot_expression_comparison(ax, df):

    y_positions = np.arange(len(GENES))

    labels = [
        f"{gene['symbol']} ({gene['gene_id']})"
        for gene in GENES
    ]

    fc1 = [
        df.loc[g["gene_id"], "GSE24177_logFC"]
        for g in GENES
    ]

    fc2 = [
        df.loc[g["gene_id"], "GSE273926_published_log2FC"]
        for g in GENES
    ]

    fdr1 = [
        df.loc[g["gene_id"], "GSE24177_FDR"]
        for g in GENES
    ]

    fdr2 = [
        df.loc[g["gene_id"], "GSE273926_published_FDR"]
        for g in GENES
    ]

    # Connect the two estimates for each gene
    for y, x1, x2 in zip(
        y_positions,
        fc1,
        fc2
    ):
        ax.plot(
            [x1, x2],
            [y, y],
            linewidth=1,
            alpha=0.3
        )

    # Larger points indicate FDR < 0.05
    size1 = [
        95 if p < 0.05 else 42
        for p in fdr1
    ]

    size2 = [
        95 if p < 0.05 else 42
        for p in fdr2
    ]

    # Study 1: microarray / limma
    ax.scatter(
        fc1,
        y_positions,
        s=size1,
        marker="o",
        label="GSE24177 — limma",
        zorder=3
    )

    # Study 2: published bulk RNA-seq DE results
    ax.scatter(
        fc2,
        y_positions,
        s=size2,
        marker="^",
        label="GSE273926 — published Table S5 U",
        zorder=3
    )

    # Add numerical effect sizes
    for y, value in zip(y_positions, fc1):

        ax.annotate(
            f"{value:+.2f}",
            (value, y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=7
        )

    for y, value in zip(y_positions, fc2):

        ax.annotate(
            f"{value:+.2f}",
            (value, y),
            xytext=(0, -15),
            textcoords="offset points",
            ha="center",
            fontsize=7
        )

    ax.axvline(
        0,
        linestyle="--",
        linewidth=1
    )

    ax.set_yticks(
        y_positions,
        labels
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "log2 Fold Change (Drought vs Control)"
    )

    ax.set_title(
        "A. Cross-study expression patterns",
        loc="left",
        fontweight="bold"
    )

    ax.grid(
        axis="x",
        alpha=0.25
    )

    ax.legend(
        loc="lower right",
        fontsize=8,
        frameon=False
    )

    ax.text(
        0.02,
        0.02,
        "Larger marker: FDR < 0.05\n"
        "Smaller marker: FDR >= 0.05\n"
        "Circle: GSE24177 | Triangle: GSE273926",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom"
    )


# ==================================================
# 5. Functional module diagram
# ==================================================

def add_module_box(
    ax,
    x,
    y,
    width,
    height,
    title,
    body
):

    text = (
        f"{title}\n"
        f"{body}"
    )

    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=9,
        linespacing=1.5,
        bbox={
            "boxstyle": "round,pad=0.6",
            "facecolor": "none",
            "linewidth": 1.1,
        }
    )


def plot_functional_modules(ax):

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.set_title(
        "B. Candidate functional modules",
        loc="left",
        fontweight="bold"
    )

    # Top row
    add_module_box(
        ax,
        0.02,
        0.59,
        0.44,
        0.27,
        "Cu/Zn-SOD machinery",
        "CSD1 ↓   CCS1 ↓   CSD2 ↓\n"
        "Concordant in both studies"
    )

    add_module_box(
        ax,
        0.54,
        0.59,
        0.44,
        0.27,
        "Fe-SOD response",
        "FSD1 ↑\n"
        "Concordant in both studies"
    )

    # Bottom row
    add_module_box(
        ax,
        0.02,
        0.24,
        0.44,
        0.27,
        "Metal handling",
        "CCH ↑   YSL2 ↑\n"
        "Concordant in both studies"
    )

    add_module_box(
        ax,
        0.54,
        0.24,
        0.44,
        0.27,
        "Copper acquisition & regulation",
        "CITF1 ↑ (regulator)\n"
        "COPT2 ↑ (Cu transporter)\n"
        "Significant only in GSE273926"
      )

    ax.text(
        0.5,
        0.10,
        "Functional grouping is descriptive.\n"
        "No causal regulatory interactions are inferred.",
        ha="center",
        va="center",
        fontsize=8,
        style="italic"
    )


# ==================================================
# 6. Main function
# ==================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    df = load_data()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 8),
        gridspec_kw={
            "width_ratios": [1.35, 1]
        }
    )

    fig.suptitle(
     "Cross-study expression patterns of candidate "
     "ROS and metal-homeostasis genes in Arabidopsis",
     fontsize=14,
     fontweight="bold"
      )

    fig.text(
        0.5,
        0.925,
        "GSE24177: early moderate drought, microarray | "
        "GSE273926: mild drought, bulk undigested-leaf RNA-seq",
        ha="center",
        fontsize=9
    )

    plot_expression_comparison(
        axes[0],
        df
    )

    plot_functional_modules(
        axes[1]
    )

    fig.text(
        0.5,
        0.025,
        "Effect sizes are shown for descriptive cross-study comparison; "
        "platforms and drought regimes differ.\n"
        "FDR was assessed independently in each study "
        "(GSE24177: limma; GSE273926: published Table S5 U). "
        "CITF1 and COPT2 were not significant in GSE24177.",
        ha="center",
        fontsize=8.5
      )

    fig.tight_layout(
        rect=[0.02, 0.07, 0.98, 0.89]
    )

    png_path = (
        OUTPUT_DIR
        / "figure1_copper_ros_drought_signature.png"
    )

    svg_path = (
        OUTPUT_DIR
        / "figure1_copper_ros_drought_signature.svg"
    )

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        svg_path,
        bbox_inches="tight"
    )

    plt.close(fig)

    print("Figure successfully created.")
    print("PNG:", png_path)
    print("SVG:", svg_path)


if __name__ == "__main__":
    main()
