import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# t-stat ranked gene
RESULTS_DIR = "/broad/mcl/members_dir/mmurali/projects/9p21/t_stat_genelist/t_stat_results"

CELL_TYPES = [
    "Cardiac_Fibroblasts_atrial",
    "Cardiac_Fibroblasts_ventricular",
    "Endothelial_Coronary_Artery",
    "Pericytes_Adipose",
    "VSMC_pulmonary_artery"
]

threshold = 2.0
fig, axes = plt.subplots(1, len(CELL_TYPES), figsize=(20, 5), sharey=True)

for i, cell_type in enumerate(CELL_TYPES):
    df = pd.read_csv(os.path.join(RESULTS_DIR, f"{cell_type}_ranked_genes.csv"))
    t_stats = df[f"{cell_type}_t_stat"].values
    abs_t_stats = np.abs(t_stats)
    significant = abs_t_stats > threshold

    ax = axes[i]
    ax.scatter(t_stats, abs_t_stats, c=significant, cmap="coolwarm", alpha=0.7)
    ax.axhline(threshold, color="black", linestyle="dashed")
    ax.set_title(cell_type.replace("_", " "), fontsize=12)
    ax.set_xlabel("t-statistic")
    if i == 0:
        ax.set_ylabel("Absolute t-statistic")

plt.tight_layout()
for ext in ["png", "svg", "pdf"]:
    plt.savefig(os.path.join(RESULTS_DIR, f"Supp_Fig2_j_Combined_Volcano_Plots.{ext}"), dpi=900)
plt.show()