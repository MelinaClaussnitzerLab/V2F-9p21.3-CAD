#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import statsmodels.api as sm

BASE_INPUT_DIR = "/broad/mcl/members_dir/mmurali/projects/9p21/t_stat_genelist"
MERGED_FILE = os.path.join(BASE_INPUT_DIR, "allCT_basal_gene_expression.csv")
OUTPUT_DIR = os.path.join(BASE_INPUT_DIR, "t_stat_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

CELL_TYPES = [
    "Cardiac_Fibroblasts_atrial",
    "Cardiac_Fibroblasts_ventricular",
    "Endothelial_Coronary_Artery",
    "Pericytes_Adipose",
    "VSMC_pulmonary_artery",
]

# Auto-detect delimiter since input files vary
with open(MERGED_FILE, "r") as f:
    delimiter = "\t" if "\t" in f.readline() else ","

df = pd.read_csv(MERGED_FILE, sep=delimiter, engine="python", on_bad_lines="skip")
df.columns = df.columns.str.strip()

# Only expression columns are needed
expr_cols = [c for c in df.columns if c.startswith("VAL_")]
if not expr_cols:
    raise ValueError("No expression columns found.")

# Normalize ENSG IDs by removing version suffix
gene_col = df.columns[0]
df[gene_col] = df[gene_col].str.split(".").str[0]
df.set_index(gene_col, inplace=True)

expr = df[expr_cols]
expr.rename(columns={c: c.replace("VAL_", "").replace(".basal", "") for c in expr_cols}, inplace=True)
expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)

for cell in CELL_TYPES:
    if cell not in expr.columns:
        continue

    # 1 for focal, -1 for others (needed for OLS contrast)
    focal = expr[[cell]]
    non_focal = expr.drop(columns=[cell])
    labels = np.concatenate([np.ones(focal.shape[1]), -np.ones(non_focal.shape[1])])
    X = labels.reshape(-1, 1)

    t_stats = []
    for gene in expr.index:
        y = expr.loc[gene].values
        model = sm.OLS(y, X)
        result = model.fit()
        t_stats.append(result.tvalues[0])

    expr[f"{cell}_t"] = t_stats
    expr[f"{cell}_abs_t"] = np.abs(t_stats)

    ranked = expr.sort_values(by=f"{cell}_abs_t", ascending=False)
    top_n = max(1, int(0.1 * len(ranked)))

    out = ranked.iloc[:top_n][[cell, f"{cell}_t", f"{cell}_abs_t"]].reset_index()
    out.rename(columns={gene_col: "GENE_ID"}, inplace=True)

    out.to_csv(os.path.join(OUTPUT_DIR, f"{cell}_ranked_genes.csv"), index=False)