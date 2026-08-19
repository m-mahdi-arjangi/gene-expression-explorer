import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# Gene Expression Dataset
# ==========================================

gene_expression = np.array([
    [12, 15, 18, 20, 17, 16],
    [10, 11, 15, 18, 14, 12],
    [22, 25, 30, 28, 24, 26],
    [14, 16, 15, 17, 18, 16],
    [30, 35, 40, 38, 36, 34]
])

print("Gene Expression Data:\n")
print(gene_expression)

# ==========================================
# Dataset Information
# ==========================================

print("\n========== Dataset Information ==========")

print("Shape:", gene_expression.shape)
print("Dimensions:", gene_expression.ndim)
print("Data Type:", gene_expression.dtype)
print("Size:", gene_expression.size)

# ==========================================
# Gene Statistics
# ==========================================

gene_mean = np.mean(gene_expression, axis=0)

print("\n========== Gene Statistics ==========")
print("Average Gene Expression:")
print(gene_mean)

# ==========================================
# Patient Statistics
# ==========================================

patient_mean = np.mean(gene_expression, axis=1)

print("\n========== Patient Statistics ==========")
print("Average Expression Per Patient:")
print(patient_mean)

highest_patient = np.argmax(patient_mean)

print(
    f"\nPatient with Highest Average Expression: "
    f"Patient {highest_patient + 1}"
)

# ==========================================
# High Expression Values
# ==========================================

high_expression = gene_expression[gene_expression > 25]

print("\n========== High Expression Values ==========")
print(high_expression)

# ==========================================
# Data Normalization
# ==========================================

normalized_expression = (
    gene_expression - np.min(gene_expression)
) / (
    np.max(gene_expression) - np.min(gene_expression)
)

print("\n========== Normalized Data ==========")
print(normalized_expression)

print("\nMinimum Value:", np.min(normalized_expression))
print("Maximum Value:", np.max(normalized_expression))

# ==========================================
# Gene Selection
# ==========================================

selected_genes = gene_mean > 20

print("\n========== Gene Selection ==========")
print("Genes with Mean Expression > 20:")
print(selected_genes)

gene_indices = np.where(gene_mean > 20)

print("\nSelected Gene Indices:")
print(gene_indices)

# ==========================================
# Visualization
# ==========================================

gene_labels = [
    "Gene1",
    "Gene2",
    "Gene3",
    "Gene4",
    "Gene5",
    "Gene6"
]

plt.figure(figsize=(8, 5))

plt.bar(gene_labels, gene_mean)

plt.title("Average Gene Expression")
plt.xlabel("Genes")
plt.ylabel("Expression Level")

plt.tight_layout()

plt.savefig("gene_expression_plot.png")

plt.show()

print("\nGraph saved as: gene_expression_plot.png")