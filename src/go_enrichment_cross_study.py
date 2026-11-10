
from pathlib import Path
import json
import ssl
import urllib.request
import urllib.error

import pandas as pd
import truststore


# --------------------------------------------------
# 1. Project paths and API configuration
# --------------------------------------------------

ROOT_DIR = Path(__file__).resolve().parents[1]

CROSS_STUDY_DIR = ROOT_DIR / "results" / "cross_study"
OUTPUT_DIR = CROSS_STUDY_DIR / "GO_enrichment"

API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

ORGANISM = "athaliana"
SIGNIFICANCE_THRESHOLD = 0.05

ssl_context = truststore.SSLContext(
    ssl.PROTOCOL_TLS_CLIENT
)


# --------------------------------------------------
# 2. Input gene sets
# --------------------------------------------------

GENE_SET_FILES = {
    "conserved_45": (
        CROSS_STUDY_DIR
        / "conserved_drought_signature_45_genes.csv"
    ),
    "shared_up_29": (
        CROSS_STUDY_DIR
        / "shared_upregulated_drought_genes_29.csv"
    ),
    "shared_down_16": (
        CROSS_STUDY_DIR
        / "shared_downregulated_drought_genes_16.csv"
    ),
}

BACKGROUND_FILE = (
    CROSS_STUDY_DIR
    / "common_tested_gene_universe.csv"
)


# --------------------------------------------------
# 3. Helper functions
# --------------------------------------------------

def load_gene_ids(file_path: Path, column: str) -> list[str]:
    """Load and clean gene identifiers from a CSV file."""

    if not file_path.exists():
        raise FileNotFoundError(
            f"Input file not found: {file_path}"
        )

    df = pd.read_csv(file_path)

    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found in {file_path.name}. "
            f"Available columns: {df.columns.tolist()}"
        )

    genes = (
        df[column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    genes = [
        gene
        for gene in genes
        if gene and gene.lower() != "nan"
    ]

    return sorted(set(genes))


def run_gprofiler(
    query_genes: list[str],
    background_genes: list[str],
) -> dict:
    """Submit a GO:BP enrichment request using verified TLS."""

    payload = {
        "organism": ORGANISM,
        "query": query_genes,
        "sources": ["GO:BP"],
        "domain_scope": "custom_annotated",
        "background": background_genes,
        "significance_threshold_method": "fdr",
        "user_threshold": SIGNIFICANCE_THRESHOLD,
        "no_iea": False,
	"no_evidences": False,
	"all_results": False,
    }

    request = urllib.request.Request(
        API_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "gene-expression-explorer/0.1",
        },
        method="POST",
    )

    try:
        with urllib.request.urlopen(
            request,
            timeout=120,
            context=ssl_context,
        ) as response:

            response_text = response.read().decode("utf-8")
            return json.loads(response_text)

    except urllib.error.HTTPError as error:
        error_body = error.read().decode(
            "utf-8",
            errors="replace"
        )

        raise RuntimeError(
            f"g:Profiler HTTP error {error.code}: "
            f"{error_body[:2000]}"
        ) from error

    except urllib.error.URLError as error:
        raise RuntimeError(
            f"Connection error: {error.reason}"
        ) from error


def save_results(
    response_data: dict,
    output_prefix: str,
) -> int:
    """Save the complete API response and enrichment table."""

    json_path = OUTPUT_DIR / f"{output_prefix}_response.json"
    csv_path = OUTPUT_DIR / f"{output_prefix}_enrichment.csv"

    with open(json_path, "w", encoding="utf-8") as file:
        json.dump(
            response_data,
            file,
            indent=4,
            ensure_ascii=False,
        )

    results = response_data.get("result", [])

    if not isinstance(results, list):
        raise ValueError(
            "Unexpected g:Profiler response: "
            "'result' is not a list."
        )

    results_df = pd.DataFrame(results)

    # Convert list/dictionary columns into JSON strings
    # so that they can be stored safely in CSV.
    for column in results_df.columns:
        results_df[column] = results_df[column].apply(
            lambda value: (
                json.dumps(value, ensure_ascii=False)
                if isinstance(value, (list, dict))
                else value
            )
        )

    results_df.to_csv(
        csv_path,
        index=False,
        encoding="utf-8-sig",
    )

    return len(results_df)


# --------------------------------------------------
# 4. Main analysis
# --------------------------------------------------

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=== Cross-study GO:BP Enrichment ===")

    # Load common background
    background_genes = load_gene_ids(
        BACKGROUND_FILE,
        "gene_id",
    )

    background_set = set(background_genes)

    print(f"\nBackground genes: {len(background_genes)}")

    summary_rows = []
    validation_rows = []

    # Process each gene set separately
    for set_name, file_path in GENE_SET_FILES.items():

        print(f"\n--- Processing: {set_name} ---")

        input_genes = load_gene_ids(
            file_path,
            "Gene",
        )

        input_set = set(input_genes)

        query_genes = sorted(
            input_set.intersection(background_set)
        )

        excluded_genes = sorted(
            input_set.difference(background_set)
        )

        print(f"Input genes: {len(input_genes)}")
        print(f"Genes in common background: {len(query_genes)}")
        print(f"Genes outside background: {len(excluded_genes)}")

        validation_rows.append({
            "gene_set": set_name,
            "input_genes": len(input_genes),
            "genes_in_background": len(query_genes),
            "genes_outside_background": len(excluded_genes),
            "excluded_gene_ids": ";".join(excluded_genes),
        })

        if len(query_genes) < 5:
            print("Skipped: fewer than five genes in background.")
            continue

        print("Submitting request to g:Profiler...")

        response_data = run_gprofiler(
            query_genes=query_genes,
            background_genes=background_genes,
        )

        output_prefix = f"{set_name}_GO_BP"

        result_count = save_results(
            response_data,
            output_prefix,
        )

        print(f"Enrichment terms returned: {result_count}")
        print(f"Saved: {output_prefix}_response.json")
        print(f"Saved: {output_prefix}_enrichment.csv")

        summary_rows.append({
            "gene_set": set_name,
            "input_genes": len(input_genes),
            "genes_in_background": len(query_genes),
            "background_size": len(background_genes),
            "returned_terms": result_count,
            "significance_threshold": SIGNIFICANCE_THRESHOLD,
            "correction_method": "fdr",
            "source": "GO:BP",
        })

    # Save gene input validation report
    validation_df = pd.DataFrame(validation_rows)

    validation_df.to_csv(
        OUTPUT_DIR / "GO_BP_gene_input_validation.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Save overall summary
    summary_df = pd.DataFrame(summary_rows)

    summary_df.to_csv(
        OUTPUT_DIR / "GO_BP_run_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n=== Analysis completed ===")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
