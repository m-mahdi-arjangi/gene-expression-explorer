# Cross-Study Drought Analysis

## Purpose

This workflow reconstructs and validates the overlap between drought-responsive genes identified in two *Arabidopsis thaliana* transcriptomics studies:

- **GSE24177:** Microarray differential expression results generated with Limma.
- **GSE273926:** Published differential expression results for the undigested (U) sample group, extracted from the article's Supplementary Table S5.

The workflow compares gene identifiers, significant gene sets, effect directions, statistical overlap, and downstream functional enrichment. It does not perform a new DESeq2 analysis of normalized counts.

This is a cross-study comparison based on the available published and reconstructed results. It is not a formal effect-size meta-analysis.

## Source Files

### GSE24177

Input:

`results/limma_early_drought_annotated.csv`

Important columns:

- `gene_id`
- `probe_id`
- `logFC`
- `adj.P.Val`
- `gene_description`

Significance threshold:

`adj.P.Val < 0.05`

### GSE273926

Expression identifier source:

`data/processed/GSE273926/GSE273926_normalized_counts_all_samples.csv`

Published differential expression reference:

`results/GSE273926/article_reference/article_table_s5_U.csv`

Important Table S5 columns:

- `Gene`
- `Log2-FC`
- `P (adj.)`

Significance threshold:

`P (adj.) < 0.05`

Table S5 is a published DEG table, not a complete table of all tested genes. Its 701 rows must not be interpreted as the full DESeq2 testing universe.

The normalized expression matrix is used as an identifier source in this workflow. It is not used to rerun differential expression analysis.

## Common Gene Identifier Universe

The script constructs the operational common identifier universe as the intersection of:

1. GSE273926 expression gene IDs.
2. All nonblank annotated gene IDs in the GSE24177 Limma results.

The historical output filename is:

`common_tested_gene_universe.csv`

Its identifier column is:

`gene_id`

### Methodological Limitation

This set is an identifier-compatible comparison universe. It is not independently verified as the complete set of genes tested by DESeq2 in GSE273926.

The published Table S5 contains significant DEG entries, not the complete tested-gene results with all adjusted p-values. Therefore, the complete GSE273926 testing universe cannot be reconstructed from Table S5 alone.

Consequently, statistical overlap results should be described as conditional on this operational common identifier universe.

A definitive common tested-gene universe requires the complete tested-gene results or an independently documented tested-gene set for both studies.

## Current Reconstructed Counts

| Quantity | Count |
|---|---:|
| Operational common identifier universe | 20,409 |
| Significant GSE24177 genes | 192 |
| Published GSE273926 U Table S5 genes | 701 |
| GSE273926 U genes inside common universe | 607 |
| GSE273926 U genes outside common universe | 94 |
| Shared significant genes | 52 |
| Same-direction shared genes | 45 |
| Opposite-direction shared genes | 7 |
| Unresolved direction among shared genes | 0 |

The 94 Table S5 genes outside the common universe are present in the GSE273926 expression matrix but lack matching annotated gene identifiers in the GSE24177 Limma results.

They are not invalid genes. They are excluded from this cross-study identifier-compatible comparison because they cannot be matched within the operational common universe.

The expression identifier audit also confirmed that `ATCG...` identifiers must be retained. Do not use a restrictive `AT[1-5]G` regular expression to construct this universe, because doing so would exclude chloroplast gene identifiers.

## Statistical Methods

Let:

- `N` = size of the operational common identifier universe.
- `K1` = significant GSE24177 genes inside the common universe.
- `K2` = significant GSE273926 U genes inside the common universe.
- `k` = observed shared significant genes.

For the current analysis:

- `N = 20,409`
- `K1 = 192`
- `K2 = 607`
- `k = 52`

The 701 published GSE273926 U genes are not all used as `K2`; only the 607 genes that fall inside the operational common identifier universe are included in the overlap calculation.

### Expected Overlap

Expected overlap is calculated as:

`K1 * K2 / N`

### Fold Enrichment

Fold enrichment is calculated as:

`observed overlap / expected overlap`

### Overlap Significance

The hypergeometric test uses the upper-tail probability for observing at least `k` overlapping genes.

The overlap Fisher exact test uses a two-sided alternative.

### Directional Concordance

Directional concordance is summarized using a 2 × 2 table:

| | GSE273926 up | GSE273926 down |
|---|---:|---:|
| GSE24177 up | up_up | up_down |
| GSE24177 down | down_up | down_down |

A two-sided Fisher exact test is applied to this directional table.

Genes with a zero logFC in either study are classified as unresolved and are not included in the 2 × 2 direction test.

The same-direction percentage is calculated using all shared genes as the denominator. Unresolved genes, if present in a future analysis, must remain explicitly reported rather than being silently assigned to either direction.

## Current Statistical Results

| Statistic | Result |
|---|---:|
| Expected overlap | 5.710422 |
| Observed overlap | 52 |
| Fold enrichment | 9.106157 |
| Hypergeometric p-value | 3.373684e-35 |
| Overlap Fisher odds ratio | 13.15861 |
| Overlap Fisher p-value | 3.373684e-35 |
| Up/up | 29 |
| Up/down | 2 |
| Down/up | 5 |
| Down/down | 16 |
| Direction Fisher odds ratio | 46.4 |
| Direction Fisher p-value | 2.261238e-07 |
| Same-direction percentage | 86.538462% |

The 52 shared genes consist of:

- 29 genes upregulated in both studies.
- 16 genes downregulated in both studies.
- 7 genes with opposite differential-expression directions.
- 0 genes with unresolved direction.

These p-values are conditional on the current operational universe definition and the published DEG selection.

They should not be interpreted as proof of biological causality, as definitive evidence of a complete shared testing universe, or as the result of a cross-study meta-analysis.

## Functional Enrichment: GO Biological Process

Downstream functional enrichment is performed on three gene sets derived from the cross-study comparison:

- **Conserved drought signature:** 45 genes with the same differential-expression direction in both studies.
- **Shared upregulated genes:** 29 genes upregulated in both studies.
- **Shared downregulated genes:** 16 genes downregulated in both studies.

The enrichment workflow uses g:Profiler for *Arabidopsis thaliana* and focuses on Gene Ontology Biological Process (GO:BP) terms.

The configured significance threshold is an adjusted p-value of 0.05, with the `custom_annotated` background setting.

The enrichment results are exploratory functional summaries of the selected gene sets. They do not establish that an enriched biological process is causally responsible for drought tolerance or that every gene assigned to a term has the same biological role.

### Gene-to-Term Membership Extraction

The enrichment response includes gene-level membership information.

The extraction workflow reconstructs the relationship between each significant GO:BP term and its intersecting query genes, preserving the adjusted p-value and available gene identifiers.

The generated membership table contains 358 rows.

Its consistency is checked against the source enrichment response. The dry-run comparison reported:

| Validation | Count |
|---|---:|
| Expected rows | 358 |
| Existing rows | 358 |
| Missing rows | 0 |
| Extra rows | 0 |
| Rows with changed adjusted p-values | 0 |

This membership table is used as the input for the gene-term network visualization.

### Figure 2: GO:BP Enrichment and Gene-Term Network

Figure 2 combines enrichment summaries and gene-to-term membership information.

The enrichment panels show:

- **Panel A:** GO:BP enrichment results for the conserved drought signature.
- **Panel B:** GO:BP enrichment results for the shared upregulated gene set.
- **Panel C:** GO:BP enrichment results for the shared downregulated gene set.

The conserved and shared downregulated panels display up to 10 selected terms each.

The shared upregulated panel displays all 5 significant terms identified for that set.

Dot size in the enrichment panels represents the number of query genes intersecting each GO term.

The network panel displays selected GO:BP terms and their associated genes:

- 12 GO:BP terms.
- 27 genes.
- 80 gene-term edges.

Network representation:

- Square nodes represent GO terms.
- GO term node size reflects the adjusted-p-value significance scale (`-log10` adjusted p-value).
- Circular nodes represent genes.
- Gene node size reflects connectivity (degree), or the number of selected terms connected to that gene.
- Orange gene nodes highlight genes with degree greater than or equal to 8.
- Edges represent gene membership in a GO term. They are not gene-gene interactions or evidence of direct molecular interactions.

The network summary identified `AT1G08830`, `AT2G28190`, and `AT4G25100` among the highly connected genes in the selected network.

High connectivity indicates repeated membership across the displayed enriched terms. It should not by itself be interpreted as proof that a gene is a master regulator or causal hub.

The current enrichment summary reports:

| Gene set | Input genes | Significant GO:BP terms | Terms displayed |
|---|---:|---:|---:|
| Conserved drought signature | 45 | 52 | 10 |
| Shared upregulated | 29 | 5 | 5 |
| Shared downregulated | 16 | 50 | 10 |

Figure files:

- `results/cross_study/GO_enrichment/figures/figure2_go_bp_enrichment_gene_term_network.png`
- `results/cross_study/GO_enrichment/figures/figure2_go_bp_enrichment_gene_term_network.svg`

The PNG is suitable for viewing and sharing. The SVG preserves vector graphics for further editing and publication-oriented figure preparation.

### Reproducibility and Validation

The membership extraction script supports a dry-run mode and an explicit write mode. The plotting workflow also has dedicated unit tests.

Run membership extraction in dry-run mode:

```powershell
python .\src\extract_go_term_membership.py
```

Write the validated membership output:

```powershell
python .\src\extract_go_term_membership.py --write
```

Generate Figure 2:

```powershell
python .\src\plot_go_bp_figure.py
```

Run the membership extraction tests:

```powershell
python -m unittest discover -s tests -p "test_extract_go_term_membership.py" -v
```

Run the Figure 2 plotting tests:

```powershell
python -m unittest discover -s tests -p "test_plot_go_bp_figure.py" -v
```

The membership extraction tests validate source parsing, gene-to-term relationships, and output consistency.

The plotting tests validate input loading and figure-generation behavior.

These tests do not replace biological validation of the enrichment findings.

## Generated Cross-Study Outputs

The primary cross-study comparison outputs are written under:

`results/cross_study/`

| File | Purpose |
|---|---|
| `common_tested_gene_universe.csv` | Operational common identifier universe |
| `cross_study_overlap_summary.csv` | Counts and directional summary |
| `cross_study_statistical_validation.csv` | Overlap and direction statistics |
| `significant_drought_genes_shared.csv` | All shared significant genes |
| `conserved_drought_signature_45_genes.csv` | Same-direction shared genes |
| `shared_upregulated_drought_genes_29.csv` | Upregulated in both studies |
| `shared_downregulated_drought_genes_16.csv` | Downregulated in both studies |
| `discordant_drought_genes_7.csv` | Opposite-direction genes |

The GO enrichment workflow also generates enrichment result tables and gene-to-term membership data under:

`results/cross_study/GO_enrichment/`

The figure files are stored in its `figures/` subdirectory.

These outputs support downstream functional interpretation and visualization.

## Running the Cross-Study Analysis

Run commands from the repository root.

### Dry Run (Default)

```powershell
python .\src\cross_study_analysis.py
```

The dry run reconstructs the outputs and compares them with the saved CSV files. It does not intentionally write or overwrite result files.

### Write Mode

```powershell
python .\src\cross_study_analysis.py --write
```

Write mode proceeds only when the existing output files match the rebuilt results. If differences are detected, it stops without overwriting the existing outputs.

Existing results should be preserved. Review any mismatch before attempting to write regenerated outputs.

## Testing

Run only the cross-study unit tests:

```powershell
python -m unittest discover -s tests -p "test_cross_study_analysis.py" -v
```

Run the complete project test suite:

```powershell
python -m unittest discover -s tests -v
```

The cross-study analysis currently has 17 dedicated unit tests.

Additional tests cover supplementary-table extraction, log2 fold-change validation, GO membership extraction, and Figure 2 generation.

The tests use small synthetic datasets where appropriate and do not depend on network access or downloading GEO files.

## Dependencies

The cross-study analysis directly uses NumPy, pandas, and SciPy.

SciPy is pinned in `requirements.txt` because the script directly uses:

- `scipy.stats.hypergeom`
- `scipy.stats.fisher_exact`

The downstream enrichment and visualization workflows have their own dependencies and external-service requirements.

Retrieving fresh g:Profiler enrichment results requires network access. The membership extraction and figure-generation steps operate on the saved enrichment outputs.

## Scientific Limitations

- GSE24177 and GSE273926 use different experimental designs and transcriptomics technologies.
- GSE273926 U results are taken from the published Table S5. Normalized counts are not used to rerun DESeq2.
- The operational common identifier universe needs further validation against the complete tested-gene sets.
- The published Table S5 contains selected DEG results, not the complete set of tested genes.
- GO enrichment depends on the selected gene sets, annotation coverage, and configured background.
- The gene-term network represents annotation membership, not physical or regulatory gene-gene interactions.
- Shared differential expression and enrichment do not establish causal mechanisms.
- The analysis is a cross-study comparison, not a formal effect-size meta-analysis.
- Additional independent datasets and sensitivity analyses are needed before making strong general claims about *Arabidopsis thaliana* drought responses.