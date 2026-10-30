# Gene Expression Explorer 🧬

**A Bioinformatics Research Project for Exploring Drought-Responsive Gene Expression in _Arabidopsis thaliana_.**

Gene Expression Explorer is a Python and R-based bioinformatics project focused on processing gene expression data, investigating biological responses to environmental stress, and identifying differentially expressed genes through statistical analysis.

The project combines data science, computational biology, and statistical bioinformatics to build a reproducible workflow for gene expression analysis.

Its long-term vision is to evolve into an extensible bioinformatics research platform integrating machine learning, functional genomics, biological pathway analysis, and explainable artificial intelligence.

---

## Research Overview

The current research focuses on understanding the early molecular response of _Arabidopsis thaliana_ to drought stress.

The project uses publicly available gene expression data from the NCBI Gene Expression Omnibus (GEO).

| Property | Description |
|---|---|
| Organism | _Arabidopsis thaliana_ |
| GEO Accession | GSE24177 |
| Platform Accession | GPL10948 |
| Research Area | Plant Molecular Biology |
| Research Focus | Drought Stress Response |
| Analysis | Differential Gene Expression |

**Dataset:** [GSE24177 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24177)

**Platform:** [GPL10948 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10948)

### Experimental Groups

The current analysis investigates the early response to moderate drought stress (Day 1).

| Experimental Group | GEO Sample IDs |
|---|---|
| Well-watered Control | GSM596123, GSM596124 |
| Moderate Drought | GSM596125, GSM596126 |

Each experimental group contains two biological replicates.

---

## Key Features

### Data Processing
- Import and process gene expression data from GEO.
- Parse GEO series matrix and platform annotation files.
- Prepare structured gene expression matrices.
- Validate expression values and sample metadata.
- Organize raw and processed datasets.

### Exploratory Data Analysis
- Analyze gene expression distributions.
- Calculate descriptive statistics.
- Investigate sample-level expression patterns.
- Perform Principal Component Analysis (PCA).
- Generate sample correlation matrices.
- Visualize gene expression profiles.

### Differential Expression Analysis
- Perform differential expression analysis using the Bioconductor `limma` package.
- Estimate gene expression differences between experimental conditions.
- Calculate moderated statistical significance.
- Apply Benjamini-Hochberg multiple testing correction.
- Identify upregulated and downregulated probes.

### Gene Annotation
- Map probe identifiers to gene identifiers.
- Retrieve gene descriptions from GEO platform annotations.
- Combine statistical results with biological annotations.
- Generate annotated candidate gene tables.

### Visualization and Reporting
- Generate expression distribution plots.
- Create expression heatmaps.
- Visualize mean expression levels.
- Generate PCA plots.
- Visualize sample correlations.
- Plot top differentially expressed candidates.
- Export analysis results as CSV files.

---

## Analysis Workflow

The project follows a reproducible multi-stage analysis workflow:

1. **Data Acquisition:** Load the original GEO dataset and platform annotation.
2. **Data Preparation:** Parse and organize the expression matrix.
3. **Data Validation:** Verify sample identifiers and expression data.
4. **Exploratory Analysis:** Investigate expression distributions and sample relationships.
5. **Dimensionality Reduction:** Apply PCA to explore overall expression patterns.
6. **Differential Expression:** Compare drought-treated samples against controls using `limma`.
7. **Multiple Testing Correction:** Calculate adjusted p-values using the Benjamini-Hochberg method.
8. **Gene Annotation:** Map probe identifiers to biological gene information.
9. **Visualization:** Generate figures and export annotated results.

---

## Current Analysis Results

The current Day 1 moderate drought comparison produced the following results:

| Analysis Metric | Result |
|---|---:|
| Total analyzed probes | 21,461 |
| Significant probes (Adjusted P-value < 0.05) | 192 |
| Upregulated probes | 108 |
| Downregulated probes | 84 |
| Probes with mapped gene identifiers | 21,389 |

The differential expression analysis uses the `limma` statistical framework.

Annotated results are available in:

- `results/limma_early_drought_results.csv`
- `results/limma_early_drought_annotated.csv`

Generated visualizations are stored in `results/figures/`.

---

## Project Structure

```text
gene-expression-explorer/
│
├── data/
│   ├── raw/
│   │   ├── GSE24177_series_matrix.txt.gz
│   │   └── GPL10948_family.soft.gz
│   │
│   ├── processed/
│   │   └── GSE24177_expression.csv
│   │
│   ├── metadata/
│   │   └── sample_metadata.csv
│   │
│   └── README.md
│
├── results/
│   ├── figures/
│   │   ├── expression_distribution.png
│   │   ├── expression_heatmap.png
│   │   ├── mean_expression.png
│   │   ├── pca_samples.png
│   │   ├── sample_correlation.png
│   │   └── top_drought_genes.png
│   │
│   ├── limma_early_drought_results.csv
│   └── limma_early_drought_annotated.csv
│
├── src/
│   ├── annotate_drought_results.py
│   ├── data_loader.py
│   ├── exploratory_analysis.py
│   ├── geo_data_loader.py
│   ├── limma_early_drought.R
│   ├── plot_drought_results.py
│   └── validate_geo_dataset.py
│
├── tests/
│   └── test_data_loader.py
│
├── install_r_dependencies.R
├── requirements.txt
├── .gitignore
└── README.md
```

---

## Technologies and Tools

| Technology | Application |
|---|---|
| Python | Data processing and analysis |
| NumPy | Numerical computation |
| Pandas | Data manipulation |
| Matplotlib | Scientific visualization |
| Scikit-learn | Principal Component Analysis |
| R | Statistical computing |
| Bioconductor limma | Differential expression analysis |
| NCBI GEO | Public gene expression datasets |
| Git & GitHub | Version control and project development |

---

## Installation

### Prerequisites

- Python 3.11 or newer
- R
- Git
- Bioconductor `limma` package

### 1. Clone the Repository

```bash
git clone https://github.com/m-mahdi-arjangi/gene-expression-explorer.git

cd gene-expression-explorer
```

### 2. Create a Python Virtual Environment

```bash
python -m venv .venv
```

Activate the environment on Windows:

```powershell
.\.venv\Scripts\Activate.ps1
```

Activate the environment on Linux or macOS:

```bash
source .venv/bin/activate
```

### 3. Install Python Dependencies

```bash
python -m pip install --upgrade pip

python -m pip install -r requirements.txt
```

The Python dependencies include:

- NumPy
- Pandas
- Matplotlib
- Scikit-learn

### 4. Install R Dependencies

Run the following command from the repository root:

```bash
Rscript install_r_dependencies.R
```

This installation script uses BiocManager to install the Bioconductor `limma` package when necessary.

---

## Running the Analysis

Run the following commands from the project root directory.

### Step 1: Prepare the GEO Dataset

```bash
python src/geo_data_loader.py
```

### Step 2: Validate the Dataset

```bash
python src/validate_geo_dataset.py
```

### Step 3: Perform Exploratory Data Analysis

```bash
python src/exploratory_analysis.py
```

### Step 4: Perform Differential Expression Analysis

```bash
Rscript src/limma_early_drought.R
```

### Step 5: Annotate Differentially Expressed Probes

```bash
python src/annotate_drought_results.py
```

### Step 6: Generate Visualizations

```bash
python src/plot_drought_results.py
```

The resulting tables and figures are saved in the `results/` directory.

---

## Running Tests

To execute the project's unit tests:

```bash
python -m unittest discover -s tests -v
```

---

## Development Roadmap

The long-term objective is to transform Gene Expression Explorer into a comprehensive bioinformatics research platform.

Planned development directions include:

### Advanced Statistical Bioinformatics
- Analysis of additional drought conditions and experimental time points.
- Advanced differential expression analysis.
- Functional enrichment analysis.
- Gene ontology and biological pathway investigation.

### Machine Learning and Artificial Intelligence
- Machine learning models for gene expression classification.
- Feature selection and candidate biomarker prioritization.
- Explainable AI for interpreting model predictions.
- Comparison of different classification algorithms.

### Systems Biology
- Gene interaction network analysis.
- Co-expression network construction.
- Identification of important regulatory genes.
- Investigation of biological pathways associated with environmental stress.

### Research and Software Development
- Validation using independent GEO datasets.
- Reproducible end-to-end analysis pipelines.
- Interactive dashboards for gene expression exploration.
- Scalable architecture for additional organisms and experiments.

---

## Scientific Reference

Harb, A., Krishnan, A., Ambavaram, M. M. R., & Pereira, A. (2010).

*Molecular and physiological analysis of drought stress in Arabidopsis reveals early responses leading to acclimation in plant growth.*

**Plant Physiology**, 154(3), 1254–1271.

https://doi.org/10.1104/pp.110.161752

Dataset: [GSE24177 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24177)

---

## Author

**Mahdi Arjangi**

Computer Engineering Student | Front-End Developer | Aspiring AI & Bioinformatics Researcher

GitHub: [@m-mahdi-arjangi](https://github.com/m-mahdi-arjangi)