# Arabidopsis Drought Transcriptomics & AI Platform

## Project Vision

Develop a reproducible bioinformatics and AI platform for studying plant stress responses using publicly available transcriptomics datasets.

The platform aims to integrate differential expression analysis, cross-study validation, functional enrichment, biological networks, machine learning, and interactive scientific visualization.

The long-term objective is to establish a reliable research foundation before introducing predictive AI models and advanced computational methods.

## Current Development Status

- [x] V1: GSE24177 drought transcriptomics analysis
- [x] V1 preservation with Git tag `v1.0.0`
- [x] V2: GSE273926 data preparation and exploratory quality control
- [x] GSE273926 Supplementary Tables S4/S5 extraction and validation
- [x] Published versus reconstructed Log2FC estimation for U, D, and FD groups
- [x] Directional agreement reconstruction with explicit zero handling
- [x] Cross-study differentially expressed gene comparison
- [x] Preliminary statistical overlap analysis
- [x] Reproducible cross-study analysis script
- [x] Cross-study output reconstruction and validation
- [x] GO:BP enrichment analysis using g:Profiler
- [x] Exploratory Figure 1: ROS and metal-homeostasis genes
- [x] GO:BP term-to-gene membership extraction and validation
- [x] Figure 2: GO functional enrichment and gene-term network
- [x] Reproducible ROS/metal candidate table generation and validation
- [x] Reproducible copper-homeostasis comparison table generation and validation
- [x] Recovery and documentation of exploratory SPL7 investigation
- [ ] Implement a standalone reproducible SPL7 analysis script
- [ ] Complete V2 repository audit and reproducibility improvements

## Development Roadmap

### Phase 1 - Repository Engineering and Reproducibility

- [ ] Audit tracked, untracked, and ignored files
- [ ] Define data storage, file-size, and version-control policies
- [ ] Review and organize source code and analysis scripts
- [x] Reproduce cross-study overlap results using a saved script
- [x] Recover GSE273926 Supplementary S4/S5 extraction workflow
- [x] Reconstruct and validate published versus estimated Log2FC outputs
- [x] Reconstruct directional agreement and explicitly handle zero-valued estimates
- [x] Recover reproducible ROS/metal and copper-homeostasis table-generation workflows
- [x] Recover and document exploratory SPL7/regulatory-gene investigation
- [ ] Implement a standalone reproducible SPL7 analysis script
- [ ] Introduce configuration files for analysis parameters and dataset settings
- [ ] Build an automated end-to-end analysis pipeline
- [ ] Add structured logging and error handling
- [ ] Document data provenance, analysis parameters, and workflow dependencies

### Phase 2 - Software Quality

- [x] Expand unit test coverage
- [ ] Add integration tests
- [x] Validate metadata and expression matrices
- [x] Test cross-study statistical calculations
- [x] Test cross-study output validation
- [x] Test Supplementary S4/S5 extraction and validation
- [x] Test Log2FC estimators, direction classification, and output comparisons
- [x] Test GO term-to-gene membership extraction and validation
- [x] Test Figure 2 generation and input validation
- [x] Test ROS/metal candidate table generation
- [x] Test copper-homeostasis table generation
- [ ] Add automated testing with GitHub Actions

**Latest verified test status: 164 unit tests passed on November 10, 2026.**

| Test Suite | Number of Tests |
|---|---:|
| Cross-study analysis | 17 |
| Data loader | 5 |
| GSE273926 preparation | 18 |
| GEO dataset validation | 8 |
| Supplementary S4/S5 extraction | 23 |
| Log2FC estimator validation | 35 |
| GO membership extraction | 12 |
| Figure 2 plotting | 8 |
| ROS/metal table generation | 15 |
| Copper-homeostasis table generation | 23 |
| **Total** | **164** |

### Phase 3 - Expanded Biological Analysis

- [ ] Integrate additional independent drought datasets
- [ ] Compare multiple tissues and drought conditions
- [ ] Perform effect-size-based meta-analysis
- [ ] Evaluate cross-study heterogeneity
- [ ] Perform sensitivity and robustness analyses
- [ ] Expand functional and regulatory network analysis
- [ ] Investigate the biological significance of conserved and discordant gene signatures

### Phase 4 - Scientific Figures and Reporting

- [x] Exploratory Figure 1: ROS and metal-homeostasis expression patterns
- [x] Figure 2: GO:BP functional enrichment visualization
- [x] Visualize genes contributing to enriched GO terms
- [ ] Prepare publication-quality figures
- [ ] Generate reproducible HTML reports
- [ ] Improve figure annotations, legends, and methodological descriptions

Figure 2 is implemented: GO:BP term-to-gene membership has been reconstructed and validated, and the enrichment visualization and gene-term network have been generated.

### Phase 5 - Machine Learning and Explainable AI

- [ ] Build a sufficiently large multi-study expression dataset
- [ ] Define appropriate prediction tasks and biological targets
- [ ] Develop drought-response classification or prediction models
- [ ] Use study-aware validation to prevent data leakage
- [ ] Evaluate models on independent datasets
- [ ] Apply feature importance and SHAP for model interpretation
- [ ] Investigate robust and reproducible gene signatures
- [ ] Compare machine learning performance against appropriate baseline models

Machine learning development should begin only after the dataset, validation strategy, and biological prediction task have been sufficiently established.

### Phase 6 - Interactive Web Platform

- [ ] Develop a Python backend using FastAPI
- [ ] Build an interactive scientific dashboard
- [ ] Develop the frontend using React and TypeScript
- [ ] Add dataset and gene search
- [ ] Visualize expression, enrichment, and network results
- [ ] Integrate validated machine learning models
- [ ] Support downloadable analysis reports
- [ ] Provide access to reproducible analysis outputs

### Phase 7 - Research Validation and Publication

- [ ] Conduct a systematic literature review
- [ ] Assess research novelty against existing studies
- [ ] Validate findings using independent datasets
- [ ] Identify opportunities for experimental validation
- [ ] Define a clear and defensible research contribution
- [ ] Prepare a scientific manuscript
- [ ] Evaluate suitable journals

## Detailed Research Execution Plan

This section expands the existing roadmap into concrete, testable research milestones. Each milestone must be completed and validated before dependent work begins.

### Research Question

Which drought-responsive gene signatures are reproducible across independent Arabidopsis transcriptomic studies, how do their expression patterns vary across tissues and cellular populations, and can these signatures support generalizable and interpretable drought-response prediction?

The project aims to investigate cross-study reproducibility, biological context, regulatory hypotheses, and predictive generalization without overstating exploratory findings.

### Milestone 1 - Complete SPL7 and Regulatory-Gene Reproducibility

- [ ] Implement a standalone SPL7 analysis script.
- [ ] Load GSE24177 annotated Limma results.
- [ ] Load GSE273926 normalized expression for the six primary undigested samples.
- [ ] Load published GSE273926 Table S5 U results.
- [ ] Reconstruct the documented SPL7 statistics and descriptive expression summaries.
- [ ] Reconstruct the curated 15-gene regulatory candidate panel.
- [ ] Validate gene identifiers, duplicate records, missing values, and required columns.
- [ ] Preserve the distinction between descriptive expression ratios and statistical log2 fold changes.
- [ ] Generate reproducible CSV outputs.
- [ ] Compare reconstructed outputs with the documented historical results.
- [ ] Add dedicated unit tests.
- [ ] Document inputs, outputs, methods, and limitations.

Completion criteria:
The historical SPL7 investigation can be reconstructed from saved inputs using a standalone script, with automated validation and no dependence on shell history.

### Milestone 2 - Configuration Management and Pipeline Automation

- [ ] Introduce centralized configuration files.
- [ ] Separate dataset paths, sample definitions, thresholds, and analysis parameters from source code.
- [ ] Define consistent input and output directory conventions.
- [ ] Build a central end-to-end analysis entry point.
- [ ] Define execution order and dependencies between analysis stages.
- [ ] Add structured logging.
- [ ] Implement clear error handling and failure reporting.
- [ ] Support dry-run and explicit write modes where appropriate.
- [ ] Protect previously validated outputs from unintended overwriting.
- [ ] Document pipeline execution and configuration options.

Completion criteria:
The major saved analysis workflows can be executed through one documented pipeline using explicit configuration and reproducible inputs.

### Milestone 3 - Independent Drought Dataset Integration

- [ ] Search GEO and other suitable public repositories for additional Arabidopsis drought transcriptomics datasets.
- [ ] Identify an initial shortlist of approximately 3-5 candidate datasets.
- [ ] Define dataset inclusion and exclusion criteria.
- [ ] Evaluate organism, tissue, developmental stage, drought treatment, duration, and biological replicates.
- [ ] Assess sequencing or microarray technology and data quality.
- [ ] Record accession numbers and original dataset provenance.
- [ ] Extract and validate sample metadata.
- [ ] Harmonize gene identifiers across studies.
- [ ] Document preprocessing and normalization procedures.
- [ ] Separate compatible and incompatible experimental comparisons.
- [ ] Preserve study-level information throughout integration.

Completion criteria:
An independently sourced multi-study dataset collection is available with validated metadata, documented preprocessing, and explicit comparability limitations.

The target number of datasets is a search objective, not a quota. Unsuitable datasets must not be included merely to increase dataset count.

### Milestone 4 - Effect-Size-Based Meta-analysis

- [ ] Define comparable biological contrasts across independent datasets.
- [ ] Obtain complete tested-gene differential expression statistics where available.
- [ ] Reconstruct or retrieve effect sizes and their standard errors.
- [ ] Validate gene-level identifier compatibility.
- [ ] Define and document the meta-analysis gene universe.
- [ ] Evaluate fixed-effect and random-effects models where assumptions permit.
- [ ] Calculate study-level heterogeneity statistics, including Q, I-squared, and tau-squared where appropriate.
- [ ] Identify genes with consistent effect directions across studies.
- [ ] Investigate genes with heterogeneous or discordant responses.
- [ ] Perform sensitivity analyses.
- [ ] Compare results under alternative inclusion and background definitions.
- [ ] Validate prioritized gene signatures using independent datasets.

Completion criteria:
A reproducible effect-size-based meta-analysis is available with documented assumptions, uncertainty estimates, heterogeneity assessment, and sensitivity analyses.

Important limitation:
Published DEG-only tables must not be treated as complete statistical result tables. Formal meta-analysis requires appropriate effect-size and uncertainty information for the genes being analyzed. The existing cross-study overlap analysis is not itself a formal meta-analysis.

### Milestone 5 - Single-cell and Cellular Context Integration

- [ ] Investigate the associated GSE273033 single-cell RNA-seq dataset.
- [ ] Review raw data availability, sample metadata, biological replicates, and experimental design.
- [ ] Establish a reproducible single-cell preprocessing workflow.
- [ ] Evaluate cell quality and filtering criteria.
- [ ] Review or validate cell-type annotations.
- [ ] Map the conserved cross-study gene signature to single-cell populations.
- [ ] Investigate cell-type-specific expression patterns.
- [ ] Compare candidate drought-responsive signatures across cellular populations.
- [ ] Evaluate pseudobulk approaches using biological replicates where appropriate.
- [ ] Compare bulk and single-cell findings.
- [ ] Investigate which expression patterns are shared or cell-type-specific.
- [ ] Document differences between bulk and single-cell measurements.

Completion criteria:
The project provides a validated cellular-context analysis of selected drought-responsive signatures, with appropriate treatment of biological replicates and cell-level observations.

### Milestone 6 - Regulatory and Network Biology

- [ ] Expand the SPL7 and regulatory-gene candidate investigation.
- [ ] Investigate HY5 and additional stress-associated transcription factors.
- [ ] Identify candidate transcription factor-target relationships using suitable external resources.
- [ ] Evaluate evidence supporting proposed regulatory relationships.
- [ ] Integrate validated gene annotations and functional relationships.
- [ ] Investigate ROS and copper-homeostasis pathways.
- [ ] Compare regulatory candidates across independent drought datasets.
- [ ] Explore co-expression and regulatory network construction where the available data support it.
- [ ] Evaluate network robustness and reproducibility.
- [ ] Distinguish experimentally supported interactions from predicted associations and correlations.
- [ ] Avoid interpreting network connectivity as proof of causality.

Completion criteria:
A documented regulatory and functional network analysis is available, with explicit evidence levels and reproducible construction procedures.

### Milestone 7 - Machine Learning and Explainable AI

Machine learning development should begin only after sufficient independent datasets, validated preprocessing, and an appropriate biological prediction task have been established.

- [ ] Define a biologically meaningful prediction target.
- [ ] Build a sufficiently large and validated multi-study expression dataset.
- [ ] Establish appropriate baseline models.
- [ ] Define study-aware data splitting strategies.
- [ ] Implement Leave-One-Study-Out validation.
- [ ] Prevent study-level data leakage during preprocessing and feature selection.
- [ ] Evaluate models on independent external datasets.
- [ ] Assess performance, calibration, and generalization.
- [ ] Investigate feature importance and SHAP.
- [ ] Evaluate stability of selected gene signatures across validation folds and studies.
- [ ] Compare predictive models with simpler statistical baselines.
- [ ] Document limitations and potential sources of bias.

Completion criteria:
Predictive models demonstrate independently evaluated performance and interpretable results without relying on random sample splitting that leaks study-specific information.

### Milestone 8 - Scientific Reporting and Research Publication

- [ ] Conduct a systematic literature review.
- [ ] Build a comparison matrix of existing studies and their contributions.
- [ ] Identify a clearly defined research gap.
- [ ] Assess scientific novelty of the proposed findings.
- [ ] Define explicit research hypotheses.
- [ ] Perform independent validation of prioritized biological findings.
- [ ] Prepare publication-quality figures.
- [ ] Generate reproducible HTML research reports.
- [ ] Document methods, datasets, statistical parameters, and limitations.
- [ ] Prepare a scientific manuscript.
- [ ] Identify appropriate journals and evaluate publication requirements.

Completion criteria:
The research contribution is clearly distinguished from prior work, independently validated where possible, and supported by reproducible analyses and transparent reporting.

### Research Dependency and Validation Gates

The following dependencies must be respected:

1. SPL7 reproducibility precedes its use in expanded regulatory analyses.
2. Configuration management and pipeline automation precede large-scale automated dataset integration.
3. Dataset quality assessment precedes cross-study integration.
4. Complete and comparable effect-size information precedes formal meta-analysis.
5. Validated cellular data and annotations precede biological interpretation of single-cell results.
6. Sufficient independent studies and a defined prediction task precede machine learning.
7. Independent validation and literature-based novelty assessment precede publication claims.

Research findings must be categorized as exploratory, computationally supported, independently validated, or experimentally validated.

No computational association, enrichment result, predictive feature, or network relationship should be presented as experimentally established causality without appropriate supporting evidence.

## Important Scientific Principles

- Public datasets must retain their original provenance.
- Exploratory estimates must not be presented as formal differential expression results.
- Published differential expression tables must not be treated as complete tested-gene result tables.
- Normalized counts must not be used as input for a new DESeq2 analysis.
- Cross-study comparisons must account for gene identifier compatibility and the definition of the tested-gene universe.
- The current operational common identifier universe must not be presented as a confirmed complete shared tested-gene universe.
- Cross-study overlap must not be confused with formal effect-size meta-analysis.
- Exact zero effect estimates must not be silently assigned to either upregulated or downregulated categories.
- Statistical significance, correlation, and functional enrichment do not establish causality.
- Machine learning models must be evaluated without study-level data leakage.
- Independent datasets should be used to evaluate the generalizability of biological findings and predictive models.
- Biological interpretations must distinguish exploratory observations from experimentally validated mechanisms.
- Research novelty must be assessed against existing scientific literature.
- All reported results and analysis workflows should be reproducible.