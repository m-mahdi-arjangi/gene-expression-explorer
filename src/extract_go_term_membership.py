from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


# --------------------------------------------------
# 1. Project configuration
# --------------------------------------------------

ROOT_DIR = Path(__file__).resolve().parents[1]

OUTPUT_DIR = ROOT_DIR / "results" / "cross_study" / "GO_enrichment"

OUTPUT_FILE = OUTPUT_DIR / "GO_BP_term_gene_membership.csv"

GENE_SETS = [
    ("Conserved 45", "conserved_45_GO_BP_response.json"),
    ("Shared up 29", "shared_up_29_GO_BP_response.json"),
    ("Shared down 16", "shared_down_16_GO_BP_response.json"),
]

OUTPUT_COLUMNS = [
    "gene_set",
    "GO_ID",
    "GO_term",
    "adjusted_p_value",
    "gene_id",
    "evidence_codes",
]


# --------------------------------------------------
# 2. Validation helpers
# --------------------------------------------------

def validate_response_settings(
    response: dict,
    filename: str,
) -> None:
    """Validate that the saved response used the expected analysis settings."""

    meta = response.get("meta", {})
    query_metadata = meta.get("query_metadata", {})

    if query_metadata.get("organism") != "athaliana":
        raise ValueError(
            f"{filename}: unexpected organism."
        )

    if query_metadata.get("significance_threshold_method") != "fdr":
        raise ValueError(
            f"{filename}: expected FDR significance correction."
        )

    threshold = query_metadata.get("user_threshold")

    if threshold is None or not math.isclose(
        float(threshold),
        0.05,
        rel_tol=0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            f"{filename}: expected significance threshold 0.05, "
            f"received {threshold}."
        )

    sources = query_metadata.get("sources", [])

    if sources != ["GO:BP"]:
        raise ValueError(
            f"{filename}: expected sources ['GO:BP'], "
            f"received {sources}."
        )


def get_query_genes(
    response: dict,
    filename: str,
) -> tuple[list[str], list[str]]:
    """
    Return the original query genes and the mapped ENSGs.

    The intersections array is indexed using the order of `ensgs`,
    not the order or length of query_size in an enrichment result.
    """

    meta = response.get("meta", {})

    query_metadata = meta.get("query_metadata", {})
    original_genes = (
        query_metadata
        .get("queries", {})
        .get("query_1")
    )

    genes_metadata = meta.get("genes_metadata", {})
    mapped_genes = (
        genes_metadata
        .get("query", {})
        .get("query_1", {})
        .get("ensgs")
    )

    if not isinstance(original_genes, list):
        raise ValueError(
            f"{filename}: query_metadata.queries.query_1 "
            "must be a list."
        )

    if not isinstance(mapped_genes, list):
        raise ValueError(
            f"{filename}: genes_metadata.query.query_1.ensgs "
            "must be a list."
        )

    if not all(
        isinstance(gene, str) and gene.strip()
        for gene in original_genes
    ):
        raise ValueError(
            f"{filename}: invalid original gene identifiers."
        )

    if not all(
        isinstance(gene, str) and gene.strip()
        for gene in mapped_genes
    ):
        raise ValueError(
            f"{filename}: invalid mapped gene identifiers."
        )

    if len(original_genes) != len(set(original_genes)):
        raise ValueError(
            f"{filename}: duplicate identifiers in original query."
        )

    if len(mapped_genes) != len(set(mapped_genes)):
        raise ValueError(
            f"{filename}: duplicate identifiers in mapped ENSGs."
        )

    unexpected = set(mapped_genes).difference(original_genes)

    if unexpected:
        raise ValueError(
            f"{filename}: mapped identifiers not found in original "
            f"query: {sorted(unexpected)[:5]}"
        )

    return original_genes, mapped_genes


# --------------------------------------------------
# 3. Membership extraction
# --------------------------------------------------

def extract_membership_rows(
    response: dict,
    gene_set: str,
    filename: str = "response",
) -> list[dict]:
    """
    Extract gene-to-GO-term membership from a g:Profiler response.

    Each non-empty intersections[i] entry belongs to ensgs[i].
    Evidence codes are preserved in their original order.
    """

    if not isinstance(response, dict):
        raise ValueError(
            f"{filename}: response must be a dictionary."
        )

    validate_response_settings(response, filename)

    original_genes, mapped_genes = get_query_genes(
        response,
        filename,
    )

    results = response.get("result")

    if not isinstance(results, list):
        raise ValueError(
            f"{filename}: result must be a list."
        )

    rows = []

    for term_index, term in enumerate(results):

        if not isinstance(term, dict):
            raise ValueError(
                f"{filename}: result[{term_index}] is not a dictionary."
            )

        if term.get("source") != "GO:BP":
            continue

        intersections = term.get("intersections")

        if not isinstance(intersections, list):
            raise ValueError(
                f"{filename}: invalid intersections for "
                f"{term.get('native')}."
            )

        if len(intersections) != len(mapped_genes):
            raise ValueError(
                f"{filename}: intersection length mismatch for "
                f"{term.get('native')}: "
                f"{len(intersections)} intersections versus "
                f"{len(mapped_genes)} mapped genes."
            )

        members = []

        for gene_id, evidence in zip(
            mapped_genes,
            intersections,
        ):

            if not isinstance(evidence, list):
                raise ValueError(
                    f"{filename}: invalid evidence list for "
                    f"{term.get('native')} / {gene_id}."
                )

            if not all(
                isinstance(code, str) and code.strip()
                for code in evidence
            ):
                raise ValueError(
                    f"{filename}: invalid evidence code for "
                    f"{term.get('native')} / {gene_id}."
                )

            if len(evidence) != len(set(evidence)):
                raise ValueError(
                    f"{filename}: duplicate evidence codes for "
                    f"{term.get('native')} / {gene_id}."
                )

            if evidence:
                members.append({
                    "gene_id": gene_id,
                    "evidence_codes": ";".join(evidence),
                })

        expected_intersection_size = term.get("intersection_size")

        if (
            not isinstance(expected_intersection_size, int)
            or isinstance(expected_intersection_size, bool)
        ):
            raise ValueError(
                f"{filename}: invalid intersection_size for "
                f"{term.get('native')}."
            )

        if len(members) != expected_intersection_size:
            raise ValueError(
                f"{filename}: intersection_size mismatch for "
                f"{term.get('native')}: "
                f"reported {expected_intersection_size}, "
                f"reconstructed {len(members)}."
            )

        # Keep membership only for significant GO:BP terms.
        if not term.get("significant", False):
            continue

        go_id = term.get("native")
        go_term = term.get("name")
        adjusted_p_value = term.get("p_value")

        if not isinstance(go_id, str) or not go_id:
            raise ValueError(
                f"{filename}: missing GO identifier."
            )

        if not isinstance(go_term, str) or not go_term:
            raise ValueError(
                f"{filename}: missing GO term name for {go_id}."
            )

        try:
            adjusted_p_value = float(adjusted_p_value)
        except (TypeError, ValueError):
            raise ValueError(
                f"{filename}: invalid p_value for {go_id}."
            )

        if (
            not math.isfinite(adjusted_p_value)
            or not 0 <= adjusted_p_value <= 1
        ):
            raise ValueError(
                f"{filename}: p_value outside [0, 1] for {go_id}."
            )

        for member in members:
            rows.append({
                "gene_set": gene_set,
                "GO_ID": go_id,
                "GO_term": go_term,
                "adjusted_p_value": adjusted_p_value,
                "gene_id": member["gene_id"],
                "evidence_codes": member["evidence_codes"],
            })

    return rows


# --------------------------------------------------
# 4. Build the complete membership table
# --------------------------------------------------

def build_membership_table(
    output_dir: Path = OUTPUT_DIR,
) -> tuple[pd.DataFrame, list[dict]]:
    """Read all saved responses and reconstruct the membership table."""

    all_rows = []
    summaries = []

    for gene_set, filename in GENE_SETS:

        json_path = output_dir / filename

        if not json_path.exists():
            raise FileNotFoundError(
                f"Required response file not found: {json_path}"
            )

        response = json.loads(
            json_path.read_text(encoding="utf-8")
        )

        original_genes, mapped_genes = get_query_genes(
            response,
            filename,
        )

        rows = extract_membership_rows(
            response=response,
            gene_set=gene_set,
            filename=filename,
        )

        all_rows.extend(rows)

        summaries.append({
            "gene_set": gene_set,
            "input_genes": len(original_genes),
            "mapped_genes": len(mapped_genes),
            "returned_terms": len(response["result"]),
            "membership_rows": len(rows),
            "unique_GO_terms": len({
                row["GO_ID"] for row in rows
            }),
        })

    dataframe = pd.DataFrame(
        all_rows,
        columns=OUTPUT_COLUMNS,
    )

    duplicate_columns = [
        "gene_set",
        "GO_ID",
        "gene_id",
    ]

    duplicates = dataframe.duplicated(
        subset=duplicate_columns,
        keep=False,
    )

    if duplicates.any():
        examples = dataframe.loc[
            duplicates,
            duplicate_columns,
        ].head(10)

        raise ValueError(
            "Duplicate gene-term memberships detected:\n"
            + examples.to_string(index=False)
        )

    return dataframe, summaries


# --------------------------------------------------
# 5. Full CSV comparison
# --------------------------------------------------

def dataframe_to_records(
    dataframe: pd.DataFrame,
) -> dict[tuple, float]:
    """Convert a membership table into comparable keyed records."""

    missing_columns = set(OUTPUT_COLUMNS).difference(
        dataframe.columns
    )

    if missing_columns:
        raise ValueError(
            f"Missing columns: {sorted(missing_columns)}"
        )

    records = {}

    for row in dataframe[OUTPUT_COLUMNS].itertuples(
        index=False,
        name=None,
    ):

        gene_set, go_id, go_term, p_value, gene_id, evidence = row

        key = (
            str(gene_set),
            str(go_id),
            str(go_term),
            str(gene_id),
            "" if pd.isna(evidence) else str(evidence),
        )

        if pd.isna(p_value):
            raise ValueError(
                f"Missing adjusted p-value for {key}."
            )

        if key in records:
            raise ValueError(
                f"Duplicate membership record: {key}"
            )

        records[key] = float(p_value)

    return records


def compare_membership_tables(
    expected: pd.DataFrame,
    existing: pd.DataFrame,
) -> dict:
    """Compare all membership records, including p-values."""

    expected_records = dataframe_to_records(expected)
    existing_records = dataframe_to_records(existing)

    expected_keys = set(expected_records)
    existing_keys = set(existing_records)

    missing = expected_keys - existing_keys
    extra = existing_keys - expected_keys

    common = expected_keys.intersection(existing_keys)

    changed_p_values = []

    for key in common:
        expected_value = expected_records[key]
        existing_value = existing_records[key]

        if not math.isclose(
            expected_value,
            existing_value,
            rel_tol=1e-12,
            abs_tol=1e-15,
        ):
            changed_p_values.append({
                "key": key,
                "expected": expected_value,
                "existing": existing_value,
            })

    return {
        "expected_rows": len(expected_records),
        "existing_rows": len(existing_records),
        "missing": missing,
        "extra": extra,
        "changed_p_values": changed_p_values,
        "matches": (
            not missing
            and not extra
            and not changed_p_values
        ),
    }


# --------------------------------------------------
# 6. Main
# --------------------------------------------------

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct GO:BP gene-term membership "
            "from saved g:Profiler JSON responses."
        )
    )

    parser.add_argument(
        "--write",
        action="store_true",
        help=(
            "Write the reconstructed membership CSV. "
            "Without this flag, run a read-only dry-run."
        ),
    )

    args = parser.parse_args()

    print("=== GO:BP Gene-Term Membership Extraction ===")
    print("Source: saved g:Profiler JSON files")
    print("Correction method: FDR")
    print("Significance threshold: 0.05")
    print("Mode:", "WRITE" if args.write else "DRY-RUN")

    expected, summaries = build_membership_table()

    print("\n--- Extraction summary ---")

    for summary in summaries:
        print(
            f"{summary['gene_set']}: "
            f"{summary['input_genes']} input genes, "
            f"{summary['mapped_genes']} mapped genes, "
            f"{summary['returned_terms']} terms, "
            f"{summary['unique_GO_terms']} significant GO terms, "
            f"{summary['membership_rows']} memberships"
        )

    print(f"\nTotal membership rows: {len(expected)}")

    if OUTPUT_FILE.exists():

        existing = pd.read_csv(
            OUTPUT_FILE,
            encoding="utf-8-sig",
        )

        comparison = compare_membership_tables(
            expected,
            existing,
        )

        print("\n--- Existing CSV comparison ---")
        print("Expected rows:", comparison["expected_rows"])
        print("Existing rows:", comparison["existing_rows"])
        print("Missing rows:", len(comparison["missing"]))
        print("Extra rows:", len(comparison["extra"]))
        print(
            "Changed p-values:",
            len(comparison["changed_p_values"]),
        )

        if comparison["matches"]:
            print("Result: FULL MATCH")
        else:
            print("Result: DIFFERENT")

            if comparison["missing"]:
                print(
                    "\nMissing examples:",
                    list(sorted(comparison["missing"]))[:5],
                )

            if comparison["extra"]:
                print(
                    "\nExtra examples:",
                    list(sorted(comparison["extra"]))[:5],
                )

            if comparison["changed_p_values"]:
                print(
                    "\nChanged p-value examples:",
                    comparison["changed_p_values"][:5],
                )

    else:
        print("\nExisting membership CSV not found.")

    if args.write:

        OUTPUT_DIR.mkdir(
            parents=True,
            exist_ok=True,
        )

        temporary_file = OUTPUT_FILE.with_suffix(
            OUTPUT_FILE.suffix + ".tmp"
        )

        expected.to_csv(
            temporary_file,
            index=False,
            encoding="utf-8-sig",
        )

        temporary_file.replace(OUTPUT_FILE)

        print("\nMembership CSV successfully written:")
        print(OUTPUT_FILE)

        print("Final rows:", len(expected))
        print(
            "Final columns:",
            expected.columns.tolist(),
        )

    else:
        print("\nDry-run completed. No files were modified.")
        print("Use --write to replace the membership CSV.")


if __name__ == "__main__":
    main()
