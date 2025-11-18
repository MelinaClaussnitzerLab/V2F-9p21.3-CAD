# MAGMA Gene-Based Analysis for GWAS

This repository contains the scripts and pipeline used to perform a MAGMA gene-based analysis on GWAS summary statistics for Coronary Artery Disease (CAD), Type 2 Diabetes (T2D), and Schizophrenia (SCZ).

The workflow is divided into two main parts:
1.  **Gene-Based Analysis**: A shell script (`01_run_magma.sh`) runs the core MAGMA analysis to generate gene-level results (`.genes.out` files) from GWAS summary statistics.
2.  **Post-processing and Plotting**: A Jupyter Notebook (`02_9p21_MAGMA_results_plotting.ipynb`) processes these output files to create publication-quality plots and supplementary tables.


## Step 1: Install Dependencies
See MAGMA github for instructions

## Step 2: Configure and Run MAGMA Analysis
```bash
bash 01_run_magma.sh
```
This will run the SNP annotation and then perform the gene-based analysis for all three traits.

## Step 3: Post-process and Plot Results
Plotting in the `02_9p21_MAGMA_results_plotting.ipynb` notebook (Fig3. D-E; Supp Fig8. A)
