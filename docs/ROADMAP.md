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