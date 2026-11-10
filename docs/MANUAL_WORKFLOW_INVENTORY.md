# Manual Workflow and Technical Debt Inventory

This document tracks analyses that were previously performed through temporary PowerShell/Python commands and records whether they have been converted into reproducible project workflows.

## A. Existing Saved Scripts and Recovered Workflows

- [x] GSE273926 data preparation
- [x] GSE273926 exploratory analysis
- [x] GSE273926 Supplementary Tables S4/S5 extraction and validation
- [x] Published versus reconstructed Log2FC validation for U, D, and FD
- [x] GO API execution
- [x] Exploratory Figure 1 generation
- [x] Cross-study overlap and statistical validation
- [x] Reproducible ROS/metal candidate table generation
- [x] Reproducible copper-homeostasis comparison tables

Core scripts recovered or independently validated during this work:

- `src/extract_gse273926_supplementary.py`
- `src/validate_log2fc_estimators.py`
- `src/cross_study_analysis.py`
- `src/extract_go_term_membership.py`
- `src/plot_go_bp_figure.py`
- `src/build_ros_metal_gene_table.py`
- `src/build_copper_homeostasis_tables.py`

## B. Manual Workflows to Recover

- [x] Extract and document published Supplementary Tables S4 and S5
- [x] Compare published and reconstructed Log2FC estimates
- [x] Validate directional agreement between estimators, including explicit zero handling
- [x] Generate cross-study DEG overlap tables
- [x] Calculate hypergeometric and Fisher statistical tests
- [x] Generate conserved, upregulated, downregulated, and discordant gene sets
- [x] Extract and independently validate GO term-to-gene membership from saved g:Profiler JSON responses
- [x] Generate GO term-to-gene membership tables
- [x] Generate Figure 2 GO enrichment visualization and gene-term network
- [x] Recover reproducible generation of ROS and copper-homeostasis candidate tables
- [x] Recover reproducible generation of copper-homeostasis sample-level comparison tables
- [x] Document exploratory SPL7 and regulatory-gene investigations

The supplementary extraction workflow reconstructs 32,833 Table S4 gene rows, validates all 14 sample mappings, checks leaf-transcriptome membership, and compares the extracted Table S5 sheets with the previously saved CSV outputs.

The Log2FC workflow reconstructs two estimators from Table S4 using a pseudocount of one. It validates the estimates against published Table S5 results, calculates correlation and error metrics, reconstructs direction comparisons, and checks existing output files before writing.

Some GO, ROS, and copper-homeostasis result tables are already present under `results/cross_study/GO_enrichment/`.

The GO membership extraction and Figure 2 visualization workflows have been recovered into independently tested scripts.

The ROS/metal candidate table workflow has been recovered into `src/build_ros_metal_gene_table.py`. Its reconstructed output was compared with the existing CSV, and all eight selected genes and nine output columns matched exactly.

The copper-homeostasis comparison workflows have been recovered into `src/build_copper_homeostasis_tables.py`. Both reconstructed outputs matched their existing CSV files exactly:

- `copper_homeostasis_regulon_comparison.csv`
- `copper_homeostasis_sample_level.csv`

The existing result files were preserved during validation.

### Recovered SPL7 and Regulatory-Gene Investigation

The exploratory SPL7 investigation focuses on Arabidopsis gene `AT5G18830` and combines three sources.

1. **GSE24177 annotated Limma results**
   - Source: `results/limma_early_drought_annotated.csv`
   - Extracted statistics: `logFC`, `P.Value`, `adj.P.Val`, and gene description.

2. **GSE273926 normalized expression**
   - Source: `data/processed/GSE273926/GSE273926_normalized_counts_undigested.csv`
   - Samples: three control and three mild-drought undigested-leaf samples.
   - Calculated values: control mean, drought mean, and descriptive drought/control ratio.

3. **Published GSE273926 Table S5 U**
   - Source: `results/GSE273926/article_reference/article_table_s5_U.csv`
   - Checked whether SPL7 appears in the published DEG list.

#### SPL7 Findings

| Analysis | Result |
| --- | ---: |
| Gene ID | `AT5G18830` |
| GSE24177 logFC | 0.29088 |
| GSE24177 adjusted P-value | 0.27888 |
| GSE24177 significance | Not significant |
| GSE273926 control mean | 1359.67 |
| GSE273926 mild-drought mean | 1401.00 |
| Descriptive drought/control ratio | 1.03040 |
| Listed in published Table S5 U | No |

The normalized expression means indicate a small descriptive increase in the mild-drought samples. This ratio is not a log2 fold-change and does not establish statistically significant differential expression.

SPL7 was not listed in the published Table S5 U DEG table. Absence from this published DEG list should not be interpreted as absence of gene expression.

#### Regulatory-Gene Candidate Panel

A curated panel of 15 candidate genes was also investigated:

SPL7, HY5, CITF1, bHLH23, COPT1, COPT2, COPT6, FRO4, FRO5, CSD1, CSD2, CCS1, FSD1, CCH, and YSL2.

Their expression patterns were compared across GSE24177 and GSE273926, including published Table S5 U statistics where available.

The candidate list was manually curated. It should not be described as the output of an automated gene-discovery algorithm.

These analyses are exploratory and descriptive. They do not establish causal regulatory relationships or experimentally validated SPL7 activity.

The SPL7 investigation has been documented from the recovered historical analysis. It does not yet have a dedicated standalone analysis script.

## C. Tests to Complete

- [x] Five data-loader unit tests
- [x] Review existing and stale test files
- [x] Add tests for GSE273926 data preparation (18 tests)
- [x] Add metadata validation tests
- [x] Add cross-study analysis tests (17 tests)
- [x] Add statistical calculation tests
- [x] Add Supplementary S4/S5 extraction tests (23 tests)
- [x] Add Log2FC estimator and output-validation tests (35 tests)
- [x] Add GO response parsing tests (12 tests)
- [x] Add Figure 2 plotting tests (8 tests)
- [x] Add cross-study output comparison tests
- [x] Add ROS/metal table reconstruction tests (15 tests)
- [x] Add copper-homeostasis table reconstruction tests (23 tests)
- [ ] Add broader integration tests

Latest verified full test run: **164 tests passed (November 10, 2026).**

Command used to run the full test suite:

```bash
python -m unittest discover -s tests -v
```

### Test Suite Breakdown

| Test Suite | Number of Tests |
| --- | ---: |
| Cross-study analysis | 17 |
| Data loader | 5 |
| GSE273926 preparation | 18 |
| GEO dataset validation | 8 |
| Supplementary S4/S5 extraction | 23 |
| Log2FC estimator validation | 35 |
| GO membership extraction | 12 |
| Figure 2 plotting | 8 |
| ROS/metal table reconstruction | 15 |
| Copper-homeostasis table reconstruction | 23 |
| **Total tests** | **164** |

The tests use synthetic fixtures where appropriate and do not require network access or downloading GEO files.

Note: `tests/test_gene_analysis.py` was renamed to `tests/test_data_loader.py` in commit `8af0403` with no content changes. The old `.pyc` file is stale cache only.

## D. Reproducibility Requirements

- [ ] Every important result must have a corresponding source script
- [ ] Each script must document its input and output files
- [ ] Analysis parameters must be recorded
- [ ] Results must be reproducible from saved inputs
- [ ] No important analysis should depend exclusively on shell history

### Recovered Workflows and Their Main Inputs

| Workflow | Source Script | Main Inputs |
| --- | --- | --- |
| Supplementary extraction | `src/extract_gse273926_supplementary.py` | Published S4/S5 supplementary workbooks, sample metadata, and membership metadata |
| Log2FC reconstruction | `src/validate_log2fc_estimators.py` | Extracted Table S4 and Table S5 CSV files |
| Cross-study analysis | `src/cross_study_analysis.py` | GSE24177 Limma results, GSE273926 expression identifiers, and published U DEG table |
| GO membership extraction | `src/extract_go_term_membership.py` | Saved g:Profiler JSON responses and enrichment result tables |
| Figure 2 generation | `src/plot_go_bp_figure.py` | GO enrichment CSV files and GO term-to-gene membership table |
| ROS/metal candidate table | `src/build_ros_metal_gene_table.py` | Conserved 45-gene cross-study signature |
| Copper regulon comparison | `src/build_copper_homeostasis_tables.py` | GSE24177 Limma results, GSE273926 normalized expression, and published Table S5 U |
| Copper sample-level comparison | `src/build_copper_homeostasis_tables.py` | GSE273926 normalized expression and published Table S5 U |
| SPL7 exploratory investigation | Documented historical workflow; no dedicated standalone script yet | GSE24177 Limma results, GSE273926 normalized expression, and published Table S5 U |

The Supplementary extraction and Log2FC scripts compare reconstructed outputs with existing saved files before writing.

The Log2FC workflow explicitly preserves zero-valued directions as unresolved and limits automatic corrections to the documented FD zero-direction case involving `AT2G16740`.

Cross-study overlap calculations have been recovered into `src/cross_study_analysis.py`. Their inputs, outputs, parameters, statistical methods, and execution commands are documented in `docs/CROSS_STUDY_ANALYSIS.md`.

The GO membership extraction and Figure 2 generation workflows have been recovered and independently tested.

The ROS/metal and copper-homeostasis table-generation workflows have also been recovered. Their reconstructed outputs were compared against the existing saved CSV files and matched exactly.

This does not mean that all reproducibility requirements have been completed. Broader integration tests, systematic configuration management, pipeline automation, and a standalone reproducible SPL7 analysis script remain outstanding.

## E. Deferred Research Deliverables

- [x] Figure 2: GO functional enrichment visualization
- [x] Functional gene-term network
- [ ] Reproducible HTML research report

**Important:** Existing result files must be preserved while recovering their generating scripts. Do not regenerate or overwrite results until the source workflow has been reviewed and reconstructed outputs have been compared with saved results.