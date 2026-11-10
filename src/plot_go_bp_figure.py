from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# --------------------------------------------------
# 1. Project configuration
# --------------------------------------------------

ROOT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = (
    ROOT_DIR
    / "results"
    / "cross_study"
    / "GO_enrichment"
)

MEMBERSHIP_FILE = DATA_DIR / "GO_BP_term_gene_membership.csv"

FIGURE_DIR = DATA_DIR / "figures"

OUTPUT_PREFIX = "figure2_go_bp_enrichment_gene_term_network"

GENE_SETS = {
    "Conserved 45": {
        "prefix": "conserved_45",
        "color": "#264653",
        "top_terms": 10,
    },
    "Shared up 29": {
        "prefix": "shared_up_29",
        "color": "#2A9D8F",
        "top_terms": 10,
    },
    "Shared down 16": {
        "prefix": "shared_down_16",
        "color": "#E76F51",
        "top_terms": 10,
    },
}

NETWORK_GENE_SET = "Conserved 45"
NETWORK_TOP_TERMS = 12

MEMBERSHIP_COLUMNS = {
    "gene_set",
    "GO_ID",
    "GO_term",
    "adjusted_p_value",
    "gene_id",
    "evidence_codes",
}

ENRICHMENT_COLUMNS = {
    "native",
    "name",
    "p_value",
    "intersection_size",
    "significant",
    "source",
}


# --------------------------------------------------
# 2. Input loading and validation
# --------------------------------------------------

def parse_boolean(series: pd.Series) -> pd.Series:
    """Convert common CSV boolean representations into booleans."""

    values = series.astype(str).str.strip().str.lower()

    return values.isin(["true", "1", "yes"])


def load_membership() -> pd.DataFrame:
    """Load and validate the gene-to-GO-term membership table."""

    if not MEMBERSHIP_FILE.exists():
        raise FileNotFoundError(
            f"Membership file not found: {MEMBERSHIP_FILE}"
        )

    membership = pd.read_csv(
        MEMBERSHIP_FILE,
        encoding="utf-8-sig",
        dtype={
            "gene_set": str,
            "GO_ID": str,
            "GO_term": str,
            "gene_id": str,
            "evidence_codes": str,
        },
    )

    missing_columns = MEMBERSHIP_COLUMNS.difference(
        membership.columns
    )

    if missing_columns:
        raise ValueError(
            "Membership CSV is missing columns: "
            f"{sorted(missing_columns)}"
        )

    if membership.empty:
        raise ValueError("Membership table is empty.")

    if membership[list(MEMBERSHIP_COLUMNS)].isna().any().any():
        raise ValueError(
            "Membership table contains missing required values."
        )

    membership["adjusted_p_value"] = pd.to_numeric(
        membership["adjusted_p_value"],
        errors="raise",
    )

    if not membership["adjusted_p_value"].between(0, 1).all():
        raise ValueError(
            "Membership table contains invalid adjusted p-values."
        )

    duplicate_mask = membership.duplicated(
        subset=["gene_set", "GO_ID", "gene_id"],
        keep=False,
    )

    if duplicate_mask.any():
        examples = membership.loc[
            duplicate_mask,
            ["gene_set", "GO_ID", "gene_id"],
        ].head(10)

        raise ValueError(
            "Duplicate gene-term memberships detected:\n"
            + examples.to_string(index=False)
        )

    return membership


def load_enrichment(
    gene_set: str,
    config: dict,
    membership: pd.DataFrame,
) -> pd.DataFrame:
    """
    Load significant GO:BP results and validate them
    against the extracted membership table.
    """

    filename = (
        f"{config['prefix']}_GO_BP_enrichment.csv"
    )

    path = DATA_DIR / filename

    if not path.exists():
        raise FileNotFoundError(
            f"Enrichment file not found: {path}"
        )

    enrichment = pd.read_csv(
        path,
        encoding="utf-8-sig",
    )

    missing_columns = ENRICHMENT_COLUMNS.difference(
        enrichment.columns
    )

    if missing_columns:
        raise ValueError(
            f"{filename} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    enrichment["p_value"] = pd.to_numeric(
        enrichment["p_value"],
        errors="raise",
    )

    enrichment["intersection_size"] = pd.to_numeric(
        enrichment["intersection_size"],
        errors="raise",
    )

    enrichment["significant"] = parse_boolean(
        enrichment["significant"]
    )

    enrichment = enrichment.loc[
        (enrichment["source"] == "GO:BP")
        & enrichment["significant"]
    ].copy()

    if enrichment.empty:
        raise ValueError(
            f"No significant GO:BP terms found for {gene_set}."
        )

    if enrichment["native"].duplicated().any():
        raise ValueError(
            f"Duplicate GO IDs in enrichment file for {gene_set}."
        )

    if not enrichment["p_value"].between(0, 1).all():
        raise ValueError(
            f"Invalid p-values in enrichment file for {gene_set}."
        )

    if (enrichment["intersection_size"] < 1).any():
        raise ValueError(
            f"Invalid intersection size in {gene_set}."
        )

    gene_edges = membership.loc[
        membership["gene_set"] == gene_set
    ].copy()

    if gene_edges.empty:
        raise ValueError(
            f"No membership records found for {gene_set}."
        )

    expected_terms = set(enrichment["native"].astype(str))
    observed_terms = set(gene_edges["GO_ID"].astype(str))

    if expected_terms != observed_terms:
        missing = expected_terms - observed_terms
        extra = observed_terms - expected_terms

        raise ValueError(
            f"Enrichment/membership GO-term mismatch for {gene_set}. "
            f"Missing memberships: {sorted(missing)[:5]}; "
            f"unexpected terms: {sorted(extra)[:5]}"
        )

    enrichment_lookup = enrichment.set_index("native")

    for go_id, group in gene_edges.groupby("GO_ID"):

        if go_id not in enrichment_lookup.index:
            raise ValueError(
                f"GO term {go_id} missing from enrichment results."
            )

        expected_size = int(
            enrichment_lookup.loc[go_id, "intersection_size"]
        )

        actual_size = group["gene_id"].nunique()

        if expected_size != actual_size:
            raise ValueError(
                f"Intersection-size mismatch for "
                f"{gene_set} / {go_id}: "
                f"enrichment={expected_size}, "
                f"membership={actual_size}."
            )

        expected_p = float(
            enrichment_lookup.loc[go_id, "p_value"]
        )

        observed_p = group["adjusted_p_value"].unique()

        if len(observed_p) != 1:
            raise ValueError(
                f"Inconsistent adjusted p-values for {go_id}."
            )

        if not np.isclose(
            expected_p,
            observed_p[0],
            rtol=1e-12,
            atol=1e-15,
        ):
            raise ValueError(
                f"Adjusted p-value mismatch for "
                f"{gene_set} / {go_id}."
            )

    return enrichment.sort_values(
        "p_value",
        ascending=True,
    ).reset_index(drop=True)


# --------------------------------------------------
# 3. Enrichment dot plots
# --------------------------------------------------

def plot_enrichment_panel(
    ax,
    enrichment: pd.DataFrame,
    gene_set: str,
    config: dict,
    panel_label: str,
) -> None:
    """Draw a ranked GO:BP enrichment dot plot."""

    top = enrichment.head(
        config["top_terms"]
    ).copy()

    top = top.sort_values(
        "p_value",
        ascending=True,
    ).reset_index(drop=True)

    top["minus_log10_p"] = -np.log10(
        top["p_value"].clip(
            lower=np.finfo(float).tiny
        )
    )

    y_positions = np.arange(len(top))

    # Dot size represents the number of intersecting genes.
    sizes = 38 + top["intersection_size"].to_numpy() * 6.5

    ax.scatter(
        top["minus_log10_p"],
        y_positions,
        s=sizes,
        color=config["color"],
        alpha=0.88,
        edgecolors="white",
        linewidths=0.8,
        zorder=3,
    )

    labels = [
        f"{row['name']}\n{row['native']}"
        for _, row in top.iterrows()
    ]

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        labels,
        fontsize=7.2,
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "-log10(FDR-adjusted p-value)",
        fontsize=8,
    )

    ax.set_title(
        f"{panel_label}  {gene_set}\n"
        f"{len(enrichment)} significant terms; "
        f"showing {len(top)}",
        loc="left",
        fontsize=10,
        fontweight="bold",
        pad=10,
    )

    max_x = max(
        float(top["minus_log10_p"].max()) * 1.18,
        1.0,
    )

    ax.set_xlim(0, max_x)

    ax.grid(
        axis="x",
        linestyle="--",
        linewidth=0.6,
        alpha=0.28,
        zorder=0,
    )

    ax.tick_params(
        axis="x",
        labelsize=7.5,
    )

    ax.tick_params(
        axis="y",
        length=0,
        pad=4,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.text(
        0.98,
        0.015,
        "Dot size = genes per term",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.5,
        color="#555555",
    )


# --------------------------------------------------
# 4. Conserved signature gene-term network
# --------------------------------------------------

def plot_gene_term_network(
    ax,
    enrichment: pd.DataFrame,
    membership: pd.DataFrame,
) -> dict:
    """
    Draw a deterministic bipartite network.

    Left side: significant GO:BP terms.
    Right side: genes connected to those terms.
    """

    selected_terms = enrichment.head(
        NETWORK_TOP_TERMS
    ).copy()

    selected_term_ids = selected_terms["native"].astype(str).tolist()

    edges = membership.loc[
        (membership["gene_set"] == NETWORK_GENE_SET)
        & (membership["GO_ID"].isin(selected_term_ids))
    ].copy()

    if edges.empty:
        raise ValueError(
            "No gene-term edges found for the selected network."
        )

    # Sort terms by adjusted p-value.
    selected_terms = selected_terms.sort_values(
        "p_value",
        ascending=True,
    ).reset_index(drop=True)

    term_ids = selected_terms["native"].astype(str).tolist()

    # Gene degree = number of selected GO terms associated with a gene.
    gene_degrees = (
        edges.groupby("gene_id")["GO_ID"]
        .nunique()
    )

    gene_order = sorted(
        gene_degrees.index.tolist(),
        key=lambda gene: (
            -int(gene_degrees.loc[gene]),
            gene,
        ),
    )

    term_y = dict(
        zip(
            term_ids,
            np.linspace(0.96, 0.04, len(term_ids)),
        )
    )

    gene_y = dict(
        zip(
            gene_order,
            np.linspace(0.97, 0.03, len(gene_order)),
        )
    )

    term_x = 0.37
    gene_x = 1.02

    # Draw edges first so that nodes and labels remain readable.
    for row in edges.itertuples(index=False):

        ax.plot(
            [term_x, gene_x],
            [
                term_y[row.GO_ID],
                gene_y[row.gene_id],
            ],
            color="#9AA5B1",
            linewidth=0.75,
            alpha=0.27,
            zorder=1,
        )

    # Term nodes.
    term_sizes = (
        95
        + (
            -np.log10(
                selected_terms["p_value"].clip(
                    lower=np.finfo(float).tiny
                )
            ).to_numpy()
            * 12
        )
    )

    ax.scatter(
        [term_x] * len(term_ids),
        [term_y[go_id] for go_id in term_ids],
        s=term_sizes,
        marker="s",
        color="#264653",
        edgecolors="white",
        linewidths=0.8,
        zorder=3,
    )

    # Gene nodes. Hubs with degree >= 8 are highlighted.
    max_degree = max(
        int(gene_degrees.max()),
        1,
    )

    gene_sizes = [
        48 + int(gene_degrees.loc[gene]) * 13
        for gene in gene_order
    ]

    gene_colors = [
        "#E76F51"
        if int(gene_degrees.loc[gene]) >= 8
        else "#457B9D"
        for gene in gene_order
    ]

    ax.scatter(
        [gene_x] * len(gene_order),
        [gene_y[gene] for gene in gene_order],
        s=gene_sizes,
        c=gene_colors,
        marker="o",
        edgecolors="white",
        linewidths=0.7,
        zorder=3,
    )

    # Term labels.
    term_lookup = selected_terms.set_index("native")

    for go_id in term_ids:

        term_name = term_lookup.loc[go_id, "name"]

        label = (
            f"{term_name}  [{go_id}]"
        )

        ax.text(
            term_x - 0.035,
            term_y[go_id],
            label,
            ha="right",
            va="center",
            fontsize=7.0,
            color="#263238",
            zorder=4,
        )

    # Gene labels. Hub degree is shown next to hub IDs.
    for gene in gene_order:

        degree = int(gene_degrees.loc[gene])

        label = (
            f"{gene}  (n={degree})"
            if degree >= 8
            else gene
        )

        ax.text(
            gene_x + 0.045,
            gene_y[gene],
            label,
            ha="left",
            va="center",
            fontsize=7.2,
            fontweight="bold" if degree >= 8 else "normal",
            color="#263238",
            zorder=4,
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#264653",
            markeredgecolor="white",
            markersize=8,
            label="GO:BP term",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#457B9D",
            markeredgecolor="white",
            markersize=7,
            label="Gene",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#E76F51",
            markeredgecolor="white",
            markersize=8,
            label="High-connectivity gene (degree ≥ 8)",
        ),
        Line2D(
            [0, 1],
            [0, 0],
            color="#9AA5B1",
            linewidth=1.2,
            label="Gene-term membership",
        ),
    ]

    ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=4,
        frameon=False,
        fontsize=7.5,
        handlelength=1.6,
        columnspacing=1.8,
    )

    ax.set_title(
        "D  Gene–term membership network: conserved drought signature",
        loc="left",
        fontsize=11,
        fontweight="bold",
        pad=15,
    )

    ax.text(
        0.5,
        1.015,
        (
            f"Top {len(term_ids)} GO:BP terms  |  "
            f"{len(gene_order)} genes  |  "
            f"{len(edges)} gene-term memberships"
        ),
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#555555",
    )

    ax.set_xlim(0, 1.58)
    ax.set_ylim(-0.025, 1.025)
    ax.axis("off")

    return {
        "selected_terms": len(term_ids),
        "genes": len(gene_order),
        "edges": len(edges),
        "hub_genes": [
            gene
            for gene in gene_order
            if int(gene_degrees.loc[gene]) >= 8
        ],
        "gene_degrees": gene_degrees.to_dict(),
    }


# --------------------------------------------------
# 5. Figure generation
# --------------------------------------------------

def create_figure(
    membership: pd.DataFrame,
    enrichment_tables: dict[str, pd.DataFrame],
) -> dict:
    """Create and save Figure 2 in PNG and SVG formats."""

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.labelcolor": "#263238",
        "text.color": "#263238",
        "xtick.color": "#455A64",
        "ytick.color": "#263238",
        "svg.fonttype": "none",
    })

    fig = plt.figure(
        figsize=(18, 15),
        facecolor="white",
    )

    grid = fig.add_gridspec(
        nrows=2,
        ncols=3,
        height_ratios=[1.22, 2.0],
        hspace=0.44,
        wspace=0.46,
    )

    panel_labels = ["A", "B", "C"]

    for index, (gene_set, config) in enumerate(
        GENE_SETS.items()
    ):

        ax = fig.add_subplot(
            grid[0, index]
        )

        plot_enrichment_panel(
            ax=ax,
            enrichment=enrichment_tables[gene_set],
            gene_set=gene_set,
            config=config,
            panel_label=panel_labels[index],
        )

    network_ax = fig.add_subplot(
        grid[1, :]
    )

    network_summary = plot_gene_term_network(
        ax=network_ax,
        enrichment=enrichment_tables[NETWORK_GENE_SET],
        membership=membership,
    )

    fig.suptitle(
        "Figure 2. GO:BP enrichment and gene-term membership "
        "across the conserved Arabidopsis drought signature",
        fontsize=14,
        fontweight="bold",
        y=0.985,
    )

    fig.text(
        0.5,
        0.015,
        (
            "Enrichment was performed with g:Profiler using the "
            "custom common tested-gene background and FDR correction "
            "(threshold = 0.05). In panels A-C, dot size represents "
            "the number of query genes intersecting each GO:BP term. "
            "In panel D, square size represents -log10(FDR-adjusted "
            "p-value), circle size represents gene degree among the "
            "selected terms, and orange circles highlight genes with "
            "degree >= 8 (a descriptive visual threshold). Edges "
            "represent gene-to-GO-term annotation membership, not "
            "gene-gene interactions. Memberships were reconstructed "
            "from saved g:Profiler intersections."
        ),
        ha="center",
        va="bottom",
        fontsize=8,
        color="#455A64",
        wrap=True,
    )

    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    png_path = FIGURE_DIR / f"{OUTPUT_PREFIX}.png"
    svg_path = FIGURE_DIR / f"{OUTPUT_PREFIX}.svg"

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    fig.savefig(
        svg_path,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    return {
        "png": png_path,
        "svg": svg_path,
        "network": network_summary,
    }


# --------------------------------------------------
# 6. Main
# --------------------------------------------------

def main() -> None:

    print("=== Figure 2: GO:BP enrichment and gene-term network ===")

    print("\nLoading membership table...")

    membership = load_membership()

    print(f"Membership records: {len(membership)}")

    enrichment_tables = {}

    print("\nLoading and validating enrichment results...")

    for gene_set, config in GENE_SETS.items():

        enrichment = load_enrichment(
            gene_set=gene_set,
            config=config,
            membership=membership,
        )

        enrichment_tables[gene_set] = enrichment

        shown_terms = min(
            len(enrichment),
            config["top_terms"],
        )

        print(
            f"{gene_set}: "
            f"{len(enrichment)} significant terms, "
            f"{shown_terms} shown in panel"
        )

    print("\nGenerating figure...")

    outputs = create_figure(
        membership=membership,
        enrichment_tables=enrichment_tables,
    )

    network = outputs["network"]

    print("\n--- Network summary ---")
    print("Selected terms:", network["selected_terms"])
    print("Genes:", network["genes"])
    print("Membership edges:", network["edges"])
    print("High-connectivity genes:", network["hub_genes"])

    print("\n--- Figure files ---")
    print("PNG:", outputs["png"])
    print("SVG:", outputs["svg"])

    print("\nFigure 2 generated successfully.")


if __name__ == "__main__":
    main()
