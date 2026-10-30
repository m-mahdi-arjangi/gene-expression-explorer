# Early Drought Differential Expression Analysis
# Dataset: GSE24177
# Organism: Arabidopsis thaliana

suppressPackageStartupMessages({
  library(limma)
})

# --------------------------------------------------
# 1. Project paths
# --------------------------------------------------

project_root <- normalizePath(
  file.path(getwd()),
  winslash = "/",
  mustWork = TRUE
)

input_file <- file.path(
  project_root,
  "data",
  "processed",
  "GSE24177_expression.csv"
)

output_dir <- file.path(
  project_root,
  "results"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

output_file <- file.path(
  output_dir,
  "limma_early_drought_results.csv"
)

# --------------------------------------------------
# 2. Load expression data
# --------------------------------------------------

cat("Loading expression data...\n")

expression_df <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("probe_id" %in% colnames(expression_df))

probe_ids <- expression_df$probe_id

expression_matrix <- as.matrix(
  expression_df[, setdiff(
    colnames(expression_df),
    "probe_id"
  ), drop = FALSE]
)

rownames(expression_matrix) <- probe_ids

storage.mode(expression_matrix) <- "numeric"

cat("Probes:", nrow(expression_matrix), "\n")
cat("Samples:", ncol(expression_matrix), "\n")

# --------------------------------------------------
# 3. Select early drought samples
# --------------------------------------------------

control_samples <- c(
  "GSM596123",
  "GSM596124"
)

drought_samples <- c(
  "GSM596125",
  "GSM596126"
)

selected_samples <- c(
  control_samples,
  drought_samples
)

if (!all(selected_samples %in% colnames(expression_matrix))) {
  stop("One or more required samples are missing.")
}

expression_matrix <- expression_matrix[
  ,
  selected_samples,
  drop = FALSE
]

# Control = reference group
# Drought = comparison group

group <- factor(
  c("control", "control", "drought", "drought"),
  levels = c("control", "drought")
)

# --------------------------------------------------
# 4. Build statistical model
# --------------------------------------------------

design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)

cat("\nDesign matrix:\n")
print(design)

# --------------------------------------------------
# 5. Differential expression analysis
# --------------------------------------------------

fit <- lmFit(
  expression_matrix,
  design
)

contrast_matrix <- makeContrasts(
  drought_vs_control = drought - control,
  levels = design
)

fit <- contrasts.fit(
  fit,
  contrast_matrix
)

fit <- eBayes(fit)

# --------------------------------------------------
# 6. Extract results
# --------------------------------------------------

results <- topTable(
  fit,
  coef = "drought_vs_control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$probe_id <- rownames(results)

results <- results[
  ,
  c(
    "probe_id",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B"
  )
]

# Exclude technical control probes from the
# biological candidate-gene results.
results <- results[
  !grepl("^AFFX-", results$probe_id),
  ,
  drop = FALSE
]

write.csv(
  results,
  output_file,
  row.names = FALSE
)

# --------------------------------------------------
# 7. Summary
# --------------------------------------------------

cat("\nAnalysis completed.\n")

cat(
  "Genes with adjusted P-value < 0.05:",
  sum(results$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Results saved to:",
  output_file,
  "\n"
)

cat("\nTop 10 results:\n")

print(
  head(
    results[
      order(results$adj.P.Val, results$P.Value),
    ],
    10
  )
)
