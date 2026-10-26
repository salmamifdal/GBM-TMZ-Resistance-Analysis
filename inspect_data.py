import pandas as pd
from pathlib import Path

file_path = (
    Path(__file__).resolve().parent
    / "data"
    / "GSE327029_mRNA-12.TPM.genes.tsv"
)

data = pd.read_csv(file_path, sep="\t")

print("Number of rows and columns:", data.shape)

print("\nColumn names:")
print(data.columns.tolist())

print("\nFirst 3 genes:")
print(data.head(3).to_string(index=False))
# Define the four groups
groups = {
    "SF126_sensitive": ["SF126-1", "SF126-2", "SF126-3"],
    "SF126_resistant": ["SF126R-1", "SF126R-2", "SF126R-3"],
    "U118MG_sensitive": ["U118MG-1", "U118MG-2", "U118MG-3"],
    "U118MG_resistant": ["U118MGR-1", "U118MGR-2", "U118MGR-3"]
}

# Calculate average expression for each group
for group_name, samples in groups.items():
    data[group_name] = data[samples].mean(axis=1)

# Look at two example genes
example_genes = data[
    data["gene_name"].isin(["TSPAN6", "DPM1"])
]

print("\nAverage gene expression:")
print(example_genes[
    ["gene_name"] + list(groups.keys())
].to_string(index=False))
import numpy as np

# Compare SF126 sensitive vs resistant
data["SF126_log2FC"] = np.log2(
    (data["SF126_resistant"] + 1) /
    (data["SF126_sensitive"] + 1)
)

# Compare U118MG sensitive vs resistant
data["U118MG_log2FC"] = np.log2(
    (data["U118MG_resistant"] + 1) /
    (data["U118MG_sensitive"] + 1)
)

# Look at our example genes
examples = data[
    data["gene_name"].isin(["TSPAN6", "DPM1"])
]

print("\nGene-expression changes:")
print(examples[
    ["gene_name", "SF126_log2FC", "U118MG_log2FC"]
].to_string(index=False))
# Check whether genes change in the same direction
same_direction = (
    (
        (data["SF126_log2FC"] > 0) &
        (data["U118MG_log2FC"] > 0)
    )
    |
    (
        (data["SF126_log2FC"] < 0) &
        (data["U118MG_log2FC"] < 0)
    )
)

# Require at least a twofold change in both models
strong_change = (
    (data["SF126_log2FC"].abs() >= 1) &
    (data["U118MG_log2FC"].abs() >= 1)
)

# Require a group mean of at least 5 TPM
# in one condition of each model
expressed = (
    data[["SF126_sensitive", "SF126_resistant"]]
    .max(axis=1).ge(5)
    &
    data[["U118MG_sensitive", "U118MG_resistant"]]
    .max(axis=1).ge(5)
)

# Select candidate genes
candidates = data[
    same_direction & strong_change & expressed
].copy()

# Rank by the smaller absolute change across models
candidates["shared_change"] = candidates[
    ["SF126_log2FC", "U118MG_log2FC"]
].abs().min(axis=1)

candidates = candidates.sort_values(
    "shared_change", ascending=False
)

print("\nNumber of candidate genes:", len(candidates))

print("\nTop 10 candidate genes:")
print(candidates[
    ["gene_name", "SF126_log2FC", "U118MG_log2FC"]
].head(10).to_string(index=False))

# Save the results
results_folder = Path(__file__).resolve().parent / "results"
results_folder.mkdir(exist_ok=True)

columns = [
    "gene_id", "gene_name",
    "SF126_sensitive", "SF126_resistant",
    "U118MG_sensitive", "U118MG_resistant",
    "SF126_log2FC", "U118MG_log2FC"
]

data[columns].to_csv(
    results_folder / "all_gene_changes.csv",
    index=False
)

candidates.to_csv(
    results_folder / "candidate_genes.csv",
    index=False
)

print("\nResults saved successfully!")
import matplotlib.pyplot as plt

# Create a folder for our graphs
figures_folder = Path(__file__).resolve().parent / "figures"
figures_folder.mkdir(exist_ok=True)

# Get the top 10 candidate genes
top10 = candidates.head(10)

# Create a scatter plot of all 44 candidates
plt.figure(figsize=(9, 6))

plt.scatter(
    candidates["SF126_log2FC"],
    candidates["U118MG_log2FC"],
    alpha=0.7
)

# Label our top 10 genes
for _, gene in top10.iterrows():
    plt.annotate(
        str(gene["gene_name"]),
        (gene["SF126_log2FC"], gene["U118MG_log2FC"]),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=8
    )

plt.axhline(0, linestyle="--")
plt.axvline(0, linestyle="--")

plt.xlabel("SF126 log2 fold change")
plt.ylabel("U118MG log2 fold change")
plt.title("Shared Gene-Expression Changes in TMZ Resistance")

plt.tight_layout()

# Save the graph
plt.savefig(figures_folder / "candidate_genes.png", dpi=200)

plt.show()

print("Graph saved successfully!")
# Count candidates with increased expression in both models
increased = (
    (candidates["SF126_log2FC"] >= 1) &
    (candidates["U118MG_log2FC"] >= 1)
).sum()

# Count candidates with decreased expression in both models
decreased = (
    (candidates["SF126_log2FC"] <= -1) &
    (candidates["U118MG_log2FC"] <= -1)
).sum()

print("\nCandidates increased in both models:", increased)
print("Candidates decreased in both models:", decreased)

# Explore genes discussed in the original study
important_genes = ["KPNA3", "MGMT", "STAT3"]

selected = data[
    data["gene_name"].isin(important_genes)
].copy()

columns_to_show = [
    "gene_name",
    "SF126_sensitive",
    "SF126_resistant",
    "U118MG_sensitive",
    "U118MG_resistant",
    "SF126_log2FC",
    "U118MG_log2FC"
]

print("\nImportant genes from the published study:")
print(selected[columns_to_show].to_string(index=False))

# Save these results for our research report
selected[columns_to_show].to_csv(
    results_folder / "important_genes.csv",
    index=False
)

print("\nImportant gene results saved!")
from itertools import combinations

# Collect all 12 sample columns
samples = [
    sample
    for group in groups.values()
    for sample in group
]

# Check for missing measurements
print("\nMissing values:",
      data[samples].isna().sum().sum())

# Check for duplicate gene IDs
print("Duplicate gene IDs:",
      data["gene_id"].duplicated().sum())

# Compare samples within each group
log_data = np.log2(data[samples] + 1)

print("\nReplicate correlations:")

for group_name, group_samples in groups.items():
    for sample1, sample2 in combinations(group_samples, 2):
        correlation = log_data[sample1].corr(
            log_data[sample2]
        )

        print(
            f"{group_name}: "
            f"{sample1} vs {sample2} = "
            f"{correlation:.3f}"
        )