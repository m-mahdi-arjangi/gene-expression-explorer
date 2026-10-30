# Gene Expression Explorer 🧬

A beginner-friendly Python project for exploring a small gene-expression matrix with **Pandas, NumPy, and Matplotlib**. It demonstrates data loading, quality checks, descriptive statistics, simple gene filtering, min–max normalization, and visualisation.

> **Important:** The included CSV is a small, synthetic teaching dataset. It is not patient data, is not from a published biological study, and must not be used to draw medical or biological conclusions.

## Features

- Load expression data from a CSV file.
- Validate the expected table structure and numeric expression values.
- Reject missing, non-numeric, duplicate, or non-finite values with clear errors.
- Calculate mean and standard deviation for each gene across samples.
- Rank genes by mean expression and filter genes using a configurable threshold.
- Apply min–max normalization to expression values.
- Save a mean-expression bar chart and an expression heatmap.
- Run unit tests for the data loader.

## Project structure

```text
gene-expression-explorer/
├── data/
│   ├── sample_gene_expression.csv
│   └── README.md
├── results/
│   └── figures/             # Generated charts are saved here
├── src/
│   ├── gene_analysis.py
│   ├── data_loader.py
│   └── visualization.py
├── tests/
│   └── test_gene_analysis.py
├── .gitignore
├── requirements.txt
└── README.md
```

## Quick start

Python 3.10+ is recommended.

1. Clone this repository and open its folder.
2. (Optional but recommended) Create and activate a virtual environment.
3. Install dependencies:

   ```bash
   python -m pip install -r requirements.txt
   ```

4. Run the analysis:

   ```bash
   python src/gene_analysis.py
   ```

5. Open the generated charts in `results/figures/`.

## Input data format

The CSV file uses one row per gene and one column per sample. The first column must be named `gene_id`; remaining columns contain numeric expression values.

```csv
gene_id,Sample_1,Sample_2,Sample_3
Gene_A,12,15,18
Gene_B,10,11,15
```

Use `data/sample_gene_expression.csv` as the example input. When using a real dataset, document its source, license, preprocessing, units, and sample/group metadata.

## Example analyses

The script reports:
- dataset dimensions;
- mean and standard deviation for each gene;
- genes above a selected mean-expression threshold;
- a min–max-normalized expression matrix;
- a mean-expression bar chart and expression heatmap.

The default threshold is an exploratory filter, **not** a statistical significance test. This project does not yet perform differential-expression analysis, correct for multiple testing, or establish biomarkers.

## Run tests

From the repository root:

```bash
python -m unittest discover -s tests -v
```

## Limitations and next steps

This repository is an educational analysis workflow, not a validated biomedical tool. The current example dataset is synthetic and small. Useful next steps include using a documented public dataset, adding sample-group metadata, implementing an appropriate differential-expression workflow, and comparing results against established tools.

## Technologies

- Python
- Pandas
- NumPy
- Matplotlib
- unittest

## Author

**Mahdi Arjangi** — Computer Engineering student and Front-End Developer exploring Python, data science, and bioinformatics.

- GitHub: [@m-mahdi-arjangi](https://github.com/m-mahdi-arjangi)
