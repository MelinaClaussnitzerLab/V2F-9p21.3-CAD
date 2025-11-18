#!/bin/bash

# This script runs a gene-based analysis using MAGMA.
# 1. Annotate SNPs to genes.
# 2. Perform gene-based tests using GWAS summary statistics for multiple traits.

# --- Configuration ---
MAGMA_EXEC="/broad/mcl/members_dir/mmurali/tools/MAGMA/magma"
PROJECT_DIR=".."
REF_DIR="${PROJECT_DIR}/data/reference"
SUMSTATS_DIR="${PROJECT_DIR}/data/sumstats"
RESULTS_DIR="${PROJECT_DIR}/results/magma_output"
REF_PANEL_PREFIX="${REF_DIR}/g1000_eur/g1000_eur"

mkdir -p ${RESULTS_DIR}/CAD ${RESULTS_DIR}/T2D ${RESULTS_DIR}/SCZ

# --- Step 1: Annotate SNPs to Genes ---
ANNOT_OUT_PREFIX="${RESULTS_DIR}/snp_annot_NCBI37"
${MAGMA_EXEC} --annotate \
    --snp-loc ${REF_PANEL_PREFIX}.bim \
    --gene-loc ${REF_DIR}/NCBI37.3.gene.loc \
    --out ${ANNOT_OUT_PREFIX}

# --- Step 2: Perform Gene-Based Analysis on GWAS Data ---

# Coronary Artery Disease (CAD)
${MAGMA_EXEC} --bfile ${REF_PANEL_PREFIX} \
    --pval ${SUMSTATS_DIR}/CardioGram_2022_CAD_GWAS_BBJ_meta.tsv use=MarkerName,P ncol=N \
    --gene-annot ${ANNOT_OUT_PREFIX}.genes.annot \
    --out ${RESULTS_DIR}/CAD/CAD_GWAS

# Type 2 Diabetes (T2D)
${MAGMA_EXEC} --bfile ${REF_PANEL_PREFIX} \
    --pval ${SUMSTATS_DIR}/DIAMANTE-EUR.sumstat.txt use=rsID,Fixed-effects_p-value N=251510 \
    --gene-annot ${ANNOT_OUT_PREFIX}.genes.annot \
    --out ${RESULTS_DIR}/T2D/T2D_GWAS

# Schizophrenia (SCZ)
${MAGMA_EXEC} --bfile ${REF_PANEL_PREFIX} \
    --pval ${SUMSTATS_DIR}/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv use=SNP,PVAL ncol=NEFF \
    --gene-annot ${ANNOT_OUT_PREFIX}.genes.annot \
    --out ${RESULTS_DIR}/SCZ/SCZ_GWAS
