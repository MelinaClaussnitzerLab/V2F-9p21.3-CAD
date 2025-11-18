# LDSC Heritability analysis scripts

**Step 1** : Run script `01_compute_tstat_genelist.py` to generate t-statistics ranking from the gene expression data.

The t-statistics ranking of genes for vascular cell wall cell-types is based on 
*Finucane et al., Nat Genet 50, 621–629 (2018)* paper

**Step 2** : Run script `02_get_degs.ipynb` to obtain top 2k differentially expressed genes (DEGs) for each cell-type

**Step 3** : Run LDSC Heritability [Terra workflow](https://app.terra.bio/#workspaces/claussnitzer-fdp/PRScs_MCL/workflows/ldsc_wf/ldsc_wf) on both the t-statics and DEG set for each cell-type

**Step 4** : Run script `03_dot_plot_DEG_hg38_basal_tstat_BC.R` to generate dot plot & corresponding results table 


