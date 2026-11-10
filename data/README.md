# Data Directory

This directory contains the source, processed, supplementary, and metadata files used by **Gene Expression Explorer**.

The project integrates two publicly available *Arabidopsis thaliana* drought-stress transcriptomics studies:

1. **GSE24177** — microarray gene expression data.
2. **GSE273926** — bulk RNA-seq expression data and published supplementary results.

The two studies have different experimental designs and data-generation workflows. Their inputs are documented separately below.

## 1. Data Sources

### GSE24177 — Microarray

- **Organism:** *Arabidopsis thaliana*
- **GEO Accession:** [GSE24177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24177)
- **Platform:** [GPL10948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10948)
- **Technology:** Microarray
- **Samples:** 12

The dataset contains control and drought-stress samples collected at different time points, including Day 1, Day 10, and the wilting stage.

### GSE273926 — Bulk RNA-seq

- **Organism:** *Arabidopsis thaliana*
- **GEO Accession:** [GSE273926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273926)
- **Technology:** Bulk RNA sequencing
- **Samples represented in the project:** 14
- **Primary analysis subset:** 6 undigested leaf samples, consisting of 3 control and 3 mild-drought samples.

The project uses published supplementary expression data and reference results from this study. The primary cross-study comparison uses the undigested leaf sample group.

**Important:** GSE273926 is a bulk RNA-seq study, not a single-cell RNA-seq dataset. The project does not currently contain the raw FASTQ reads for this study.

## 2. Directory Structure

The current data directory contains:

```text
data/
├── README.md
│
├── raw/
│   ├── GSE24177_series_matrix.txt.gz
│   └── GPL10948_family.soft.gz
│
├── processed/
│   ├── GSE24177_expression.csv
│   └── GSE273926/
│       ├── GSE273926_normalized_counts_all_samples.csv
│       └── GSE273926_normalized_counts_undigested.csv
│
├── metadata/
│   ├── GSE273926_all_RunInfo.csv
│   ├── GSE273926_GEO_SRA_mapping.csv
│   ├── GSE273926_leaf_transcriptome_membership.csv
│   ├── GSE273926_sample_metadata.csv
│   ├── GSE273926_selected_RunInfo.csv
│   ├── GSE273926_selected_undigested_RunInfo.csv
│   ├── GSE273926_sra_records.xml
│   ├── GSE273926_sra_search.json
│   ├── PMC11982798_s3_listing.xml
│   └── sample_metadata.csv
│
└── supplementary/
    └── GSE273926/
        ├── GSE273926_TenorioBerrio_2024_processed.xlsx
        ├── NPH-246-840-s001.xlsx
        └── NPH-246-840-s010.xlsx
```

## 3. Raw GEO Source Files

The `raw/` directory contains the original compressed GEO source files currently used for GSE24177:

- `GSE24177_series_matrix.txt.gz` — GEO series matrix containing microarray expression measurements and sample-related information.
- `GPL10948_family.soft.gz` — GEO platform annotation containing probe and gene information.

These files are source records obtained from GEO. The series matrix contains processed microarray measurements; it should not be interpreted as raw instrument-level microarray data.

The GSE273926 workflow does not currently use locally stored raw FASTQ files.

## 4. Processed Expression Data

### GSE24177

File:

```text
data/processed/GSE24177_expression.csv
```

This expression matrix is generated from the GSE24177 GEO series matrix.

- Rows represent microarray probe identifiers.
- Columns represent GEO sample identifiers.
- Values represent the expression measurements provided by the GEO series matrix.

The matrix is used for exploratory analysis, quality assessment, and downstream differential-expression analysis.

### GSE273926

Two processed expression matrices are available.

**All samples:**

```text
data/processed/GSE273926/GSE273926_normalized_counts_all_samples.csv
```

This file contains normalized expression data for 32,833 genes across 14 samples.

**Undigested leaf subset:**

```text
data/processed/GSE273926/GSE273926_normalized_counts_undigested.csv
```

This file contains the six undigested leaf samples selected for the primary analysis: three controls and three mild-drought samples.

Both matrices contain gene-level expression data. Rows represent gene identifiers, and columns represent sample identifiers.

**Analysis limitation:** These expression matrices are already normalized. The project does not treat them as raw count data and does not rerun DESeq2 on the normalized matrices.

## 5. Sample and Study Metadata

The `metadata/` directory contains supporting files used to identify samples, track GEO/SRA relationships, and document study membership.

The files include:

| File | Role |
|---|---|
| `sample_metadata.csv` | GSE24177 sample metadata |
| `GSE273926_sample_metadata.csv` | GSE273926 sample metadata |
| `GSE273926_all_RunInfo.csv` | SRA RunInfo records associated with the study |
| `GSE273926_selected_RunInfo.csv` | Selected SRA RunInfo records |
| `GSE273926_selected_undigested_RunInfo.csv` | SRA RunInfo records for the selected undigested samples |
| `GSE273926_GEO_SRA_mapping.csv` | GEO-to-SRA sample and run mapping |
| `GSE273926_leaf_transcriptome_membership.csv` | Leaf transcriptome membership information |
| `GSE273926_sra_records.xml` | Retrieved SRA record information |
| `GSE273926_sra_search.json` | SRA search results and associated metadata |
| `PMC11982798_s3_listing.xml` | Supplementary-file listing metadata associated with the publication |

These files support sample identification, experimental-group assignment, provenance tracking, and reproducibility.

For GSE273926, the distinction between all study samples and the six undigested leaf samples is important. The primary analysis uses the latter subset.

## 6. Supplementary Data

The `supplementary/GSE273926/` directory contains published supplementary workbooks associated with the GSE273926 study:

- `GSE273926_TenorioBerrio_2024_processed.xlsx`
- `NPH-246-840-s001.xlsx`
- `NPH-246-840-s010.xlsx`

These workbooks provide supplementary expression data and supporting study information used by the GSE273926 preparation and validation workflows.

Extracted and prepared reference tables are stored separately under:

```text
results/GSE273926/article_reference/
```

The published Table S5 results are used as reference differential-expression results. They are not equivalent to a complete matrix of all genes tested in the original statistical analysis.

## 7. Data Processing and Validation

### GSE24177

The GEO expression matrix is prepared using:

```bash
python src/geo_data_loader.py
```

Dataset consistency and sample metadata are checked using:

```bash
python src/validate_geo_dataset.py
```

The resulting expression matrix is used in downstream exploratory analysis and differential-expression analysis.

### GSE273926

The GSE273926 preparation workflow uses the available published supplementary expression data and study metadata.

The workflow prepares the expression matrices for all available samples and the selected undigested subset. Supporting metadata and quality-control outputs are maintained separately.

The project preserves the distinction between:

- All samples represented in the study metadata.
- The six undigested leaf samples selected for the primary analysis.
- Published supplementary differential-expression results.

The published reference results are used as such; they are not presented as results independently recalculated from raw sequencing reads.

## 8. Data Provenance and Reproducibility

The project uses publicly available data from NCBI GEO, SRA-related metadata records, and published supplementary materials associated with GSE273926.

The provenance of each analysis input is maintained through source files, metadata, preparation scripts, and generated results.

When interpreting the outputs, the following distinctions should be maintained:

- GSE24177 provides GEO microarray series-matrix measurements and platform annotations.
- GSE273926 provides published supplementary expression data and reference results used by the current workflow.
- The primary GSE273926 cross-study comparison uses the six undigested leaf samples.
- Published differential-expression tables are reference results, not a substitute for raw sequencing data or the complete tested-gene universe.
- Normalized expression data must not be treated as raw counts for a new DESeq2 analysis.

These distinctions are important for reproducibility and for avoiding unsupported conclusions when comparing results across studies.