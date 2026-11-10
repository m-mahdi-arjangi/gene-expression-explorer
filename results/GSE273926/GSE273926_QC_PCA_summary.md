# GSE273926 – Exploratory QC and PCA Summary

## Dataset overview

- Organism: Arabidopsis thaliana
- Dataset: GSE273926
- Data type: Bulk RNA-seq
- Primary comparison: Undigested leaf samples
- Total genes: 32,833
- Total samples: 6
- Conditions: Control (3), Mild drought (3)

## Data processing

The source workbook provides normalized expression values.

For exploratory visualization, expression values were transformed using:

`log2(normalized_expression + 1)`

These values are not raw counts and were not used for DESeq2.

## Sample quality assessment

- Missing expression values: 0
- Negative expression values: 0
- Duplicate gene identifiers: 0
- Pairwise Pearson correlation: approximately 0.993–0.996

## Experimental rounds

The six primary samples belong to two experimental rounds:

- Experiment 1: one control and one mild drought sample.
- Experiment 2: two control and two mild drought samples.

The `experiment_round` variable is included in the sample metadata.

## PCA methodology

- Transformation: log2(x + 1)
- Features: 1,000 most variable genes
- Components: 2

### Explained variance

| Component | Explained variance |
|---|---:|
| PC1 | 35.01% |
| PC2 | 24.37% |
| Combined | 59.38% |

### Exploratory observations

1. Control samples have negative PC1 coordinates.
2. Mild drought samples have positive PC1 coordinates.
3. Experiment 1 samples have positive PC2 coordinates.
4. Experiment 2 samples have negative PC2 coordinates.

These patterns suggest that both experimental condition and experimental round are associated with variation in the expression matrix.

## Limitations

- Only three biological samples are available per condition.
- Experiment 1 contains one sample per condition.
- PCA is descriptive and does not establish statistical significance.
- The expression matrix contains normalized values, not raw counts.
- Differential expression analysis using DESeq2 has not been performed.

## Generated outputs

- `sample_qc_summary.csv`
- `pca_coordinates.csv`
- `pca_group_summary.csv`
- `figures/expression_distribution.png`
- `figures/sample_correlation.png`
- `figures/expression_heatmap.png`
- `figures/pca_samples_condition_and_round.png`

## Next steps

1. Investigate the contribution of experimental round to expression variation.
2. Identify an appropriate statistical analysis for normalized expression data.
3. Compare drought-associated gene signatures with GSE24177.
4. Perform pathway enrichment analysis after obtaining a defensible gene-level result.
