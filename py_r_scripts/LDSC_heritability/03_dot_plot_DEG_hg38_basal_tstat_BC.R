# --------------------------------------------------------------------
#  Dot plot code for LDSC Heritability results-  Bonferroni Corrected
# --------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(stringr)


# File paths for each cell type
cfa_paths <- c(
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/1707/Cardiac_fibroblasts_atrial_basal_vs_hg_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/T1D/Cardiac_fibroblasts_atrial_basal_vs_hg_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/aging/Cardiac_fibroblasts_atrial_basal_vs_hg_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_atrial_basal_vs_hg_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/glaucoma/Cardiac_fibroblasts_atrial_basal_vs_hg_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/ischemic_stroke/Cardiac_fibroblasts_atrial_basal_vs_hg_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/t2d_2024/Cardiac_fibroblasts_atrial_basal_vs_hg_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/t2d_hm3/Cardiac_fibroblasts_atrial_basal_vs_hg_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/total_cholesterol/Cardiac_fibroblasts_atrial_basal_vs_hg_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/triacylglycerol/Cardiac_fibroblasts_atrial_basal_vs_hg_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/1707/Cardiac_fibroblasts_atrial_basal_vs_IL1a_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/T1D/Cardiac_fibroblasts_atrial_basal_vs_IL1a_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/aging/Cardiac_fibroblasts_atrial_basal_vs_IL1a_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_atrial_basal_vs_IL1a_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/glaucoma/Cardiac_fibroblasts_atrial_basal_vs_IL1a_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/ischemic_stroke/Cardiac_fibroblasts_atrial_basal_vs_IL1a_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/t2d_2024/Cardiac_fibroblasts_atrial_basal_vs_IL1a_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/t2d_hm3/Cardiac_fibroblasts_atrial_basal_vs_IL1a_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/total_cholesterol/Cardiac_fibroblasts_atrial_basal_vs_IL1a_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/triacylglycerol/Cardiac_fibroblasts_atrial_basal_vs_IL1a_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/1707/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/T1D/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/aging/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/glaucoma/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/ischemic_stroke/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/t2d_2024/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/t2d_hm3/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/total_cholesterol/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/triacylglycerol/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/1707/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/T1D/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/aging/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/glaucoma/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/ischemic_stroke/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/t2d_2024/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/t2d_hm3/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/total_cholesterol/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/triacylglycerol/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/T1D/Cardiac_fibroblasts_atrial_basal_tstat_T1D.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/1707/Cardiac_fibroblasts_atrial_basal_tstat_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/aging/Cardiac_fibroblasts_atrial_basal_tstat_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_atrial_basal_tstat_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/glaucoma/Cardiac_fibroblasts_atrial_basal_tstat_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/ischemic_stroke/Cardiac_fibroblasts_atrial_basal_tstat_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/t2d_2024/Cardiac_fibroblasts_atrial_basal_tstat_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/t2d_hm3/Cardiac_fibroblasts_atrial_basal_tstat_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/total_cholesterol/Cardiac_fibroblasts_atrial_basal_tstat_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/triacylglycerol/Cardiac_fibroblasts_atrial_basal_tstat_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/6142/Cardiac_fibroblasts_atrial_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/Height/Cardiac_fibroblasts_atrial_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_tstat_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/6142/Cardiac_fibroblasts_atrial_basal_vs_hg_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/6142/Cardiac_fibroblasts_atrial_basal_vs_IL1a_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/6142/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/6142/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/6142/Cardiac_fibroblasts_atrial_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/Height/Cardiac_fibroblasts_atrial_basal_vs_hg_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/Height/Cardiac_fibroblasts_atrial_basal_vs_IL1a_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/Height/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/Height/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/Height/Cardiac_fibroblasts_atrial_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_hg/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_vs_hg_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL1a/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_vs_IL1a_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_IL6TNFa/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_vs_IL6TNFa_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA/basal_vs_oxLDL/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_vs_oxLDL_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFA/basal_tstat/03_h2/SCZ/Cardiac_fibroblasts_atrial_basal_tstat_SCZ.results'
)

cfv_paths <- c(
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_vs_hg_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_vs_hg_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_vs_hg_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_vs_hg_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_vs_hg_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_vs_hg_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_vs_hg_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_vs_hg_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_vs_hg_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_vs_hg_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_triacylglycerol.results',  
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/1707/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/T1D/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/aging/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/1707/Cardiac_fibroblasts_ventricular_basal_tstat_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/T1D/Cardiac_fibroblasts_ventricular_basal_tstat_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/aging/Cardiac_fibroblasts_ventricular_basal_tstat_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/CAD_GWAS_BBJ_meta_hm3/Cardiac_fibroblasts_ventricular_basal_tstat_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/glaucoma/Cardiac_fibroblasts_ventricular_basal_tstat_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/ischemic_stroke/Cardiac_fibroblasts_ventricular_basal_tstat_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/t2d_2024/Cardiac_fibroblasts_ventricular_basal_tstat_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/t2d_hm3/Cardiac_fibroblasts_ventricular_basal_tstat_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/total_cholesterol/Cardiac_fibroblasts_ventricular_basal_tstat_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/triacylglycerol/Cardiac_fibroblasts_ventricular_basal_tstat_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/CFV/basal_tstat/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_tstat_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_vs_hg_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/6142/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/6142/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_vs_hg_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/Height/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_low_serum/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_vs_low_serum_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_hg/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_vs_hg_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL1a/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_vs_IL1a_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_IL6TNFa/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_vs_IL6TNFa_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/basal_vs_oxLDL/03_h2/SCZ/Cardiac_fibroblasts_ventricular_basal_vs_oxLDL_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/Height/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFV/low_serum_vs_oxLDL/03_h2/SCZ/Cardiac_fibroblasts_ventricular_low_serum_vs_oxLDL_SCZ.results'
)

ec_paths <- c(
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/1707/Endothelial_coronary_artery_basal_vs_hg_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/T1D/Endothelial_coronary_artery_basal_vs_hg_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/aging/Endothelial_coronary_artery_basal_vs_hg_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_basal_vs_hg_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/glaucoma/Endothelial_coronary_artery_basal_vs_hg_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/ischemic_stroke/Endothelial_coronary_artery_basal_vs_hg_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/t2d_2024/Endothelial_coronary_artery_basal_vs_hg_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/t2d_hm3/Endothelial_coronary_artery_basal_vs_hg_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/total_cholesterol/Endothelial_coronary_artery_basal_vs_hg_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/triacylglycerol/Endothelial_coronary_artery_basal_vs_hg_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/1707/Endothelial_coronary_artery_basal_vs_IL6TNFa_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/T1D/Endothelial_coronary_artery_basal_vs_IL6TNFa_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/aging/Endothelial_coronary_artery_basal_vs_IL6TNFa_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_basal_vs_IL6TNFa_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/glaucoma/Endothelial_coronary_artery_basal_vs_IL6TNFa_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/ischemic_stroke/Endothelial_coronary_artery_basal_vs_IL6TNFa_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/t2d_2024/Endothelial_coronary_artery_basal_vs_IL6TNFa_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/t2d_hm3/Endothelial_coronary_artery_basal_vs_IL6TNFa_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/total_cholesterol/Endothelial_coronary_artery_basal_vs_IL6TNFa_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/triacylglycerol/Endothelial_coronary_artery_basal_vs_IL6TNFa_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/1707/Endothelial_coronary_artery_basal_vs_low_serum_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/T1D/Endothelial_coronary_artery_basal_vs_low_serum_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/aging/Endothelial_coronary_artery_basal_vs_low_serum_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_basal_vs_low_serum_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/glaucoma/Endothelial_coronary_artery_basal_vs_low_serum_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/ischemic_stroke/Endothelial_coronary_artery_basal_vs_low_serum_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/t2d_2024/Endothelial_coronary_artery_basal_vs_low_serum_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/t2d_hm3/Endothelial_coronary_artery_basal_vs_low_serum_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/total_cholesterol/Endothelial_coronary_artery_basal_vs_low_serum_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/triacylglycerol/Endothelial_coronary_artery_basal_vs_low_serum_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/1707/Endothelial_coronary_artery_basal_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/T1D/Endothelial_coronary_artery_basal_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/aging/Endothelial_coronary_artery_basal_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_basal_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/glaucoma/Endothelial_coronary_artery_basal_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/ischemic_stroke/Endothelial_coronary_artery_basal_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/t2d_2024/Endothelial_coronary_artery_basal_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/t2d_hm3/Endothelial_coronary_artery_basal_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/total_cholesterol/Endothelial_coronary_artery_basal_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/triacylglycerol/Endothelial_coronary_artery_basal_vs_oxLDL_triacylglycerol.results',  
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/1707/Endothelial_coronary_artery_low_serum_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/T1D/Endothelial_coronary_artery_low_serum_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/aging/Endothelial_coronary_artery_low_serum_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_low_serum_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/glaucoma/Endothelial_coronary_artery_low_serum_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/ischemic_stroke/Endothelial_coronary_artery_low_serum_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/t2d_2024/Endothelial_coronary_artery_low_serum_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/t2d_hm3/Endothelial_coronary_artery_low_serum_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/total_cholesterol/Endothelial_coronary_artery_low_serum_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/triacylglycerol/Endothelial_coronary_artery_low_serum_vs_oxLDL_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/1707/Endothelial_coronary_artery_basal_tstat_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/T1D/Endothelial_coronary_artery_basal_tstat_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/aging/Endothelial_coronary_artery_basal_tstat_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/CAD_GWAS_BBJ_meta_hm3/Endothelial_coronary_artery_basal_tstat_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/glaucoma/Endothelial_coronary_artery_basal_tstat_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/ischemic_stroke/Endothelial_coronary_artery_basal_tstat_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/t2d_2024/Endothelial_coronary_artery_basal_tstat_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/t2d_hm3/Endothelial_coronary_artery_basal_tstat_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/total_cholesterol/Endothelial_coronary_artery_basal_tstat_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/triacylglycerol/Endothelial_coronary_artery_basal_tstat_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/6142/Endothelial_coronary_artery_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/Height/Endothelial_coronary_artery_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/EC/basal_tstat/03_h2/SCZ/Endothelial_coronary_artery_basal_tstat_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/6142/Endothelial_coronary_artery_basal_vs_hg_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/6142/Endothelial_coronary_artery_basal_vs_IL6TNFa_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/6142/Endothelial_coronary_artery_basal_vs_low_serum_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/6142/Endothelial_coronary_artery_basal_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/6142/Endothelial_coronary_artery_low_serum_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/Height/Endothelial_coronary_artery_basal_vs_hg_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/Height/Endothelial_coronary_artery_basal_vs_IL6TNFa_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/Height/Endothelial_coronary_artery_basal_vs_oxLDL_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/Height/Endothelial_coronary_artery_basal_vs_low_serum_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_low_serum/03_h2/SCZ/Endothelial_coronary_artery_basal_vs_low_serum_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_hg/03_h2/SCZ/Endothelial_coronary_artery_basal_vs_hg_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_IL6TNFa/03_h2/SCZ/Endothelial_coronary_artery_basal_vs_IL6TNFa_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/basal_vs_oxLDL/03_h2/SCZ/Endothelial_coronary_artery_basal_vs_oxLDL_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/Height/Endothelial_coronary_artery_low_serum_vs_oxLDL_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/EC/low_serum_vs_oxLDL/03_h2/SCZ/Endothelial_coronary_artery_low_serum_vs_oxLDL_SCZ.results'
  
)

pa_paths <- c(
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/1707/Pericytes_Adipose_basal_vs_hg_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/T1D/Pericytes_Adipose_basal_vs_hg_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/aging/Pericytes_Adipose_basal_vs_hg_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_vs_hg_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/glaucoma/Pericytes_Adipose_basal_vs_hg_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/ischemic_stroke/Pericytes_Adipose_basal_vs_hg_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/t2d_2024/Pericytes_Adipose_basal_vs_hg_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/t2d_hm3/Pericytes_Adipose_basal_vs_hg_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/total_cholesterol/Pericytes_Adipose_basal_vs_hg_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/triacylglycerol/Pericytes_Adipose_basal_vs_hg_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/1707/Pericytes_Adipose_basal_vs_IL1a_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/T1D/Pericytes_Adipose_basal_vs_IL1a_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/aging/Pericytes_Adipose_basal_vs_IL1a_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_vs_IL1a_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/glaucoma/Pericytes_Adipose_basal_vs_IL1a_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/ischemic_stroke/Pericytes_Adipose_basal_vs_IL1a_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/t2d_2024/Pericytes_Adipose_basal_vs_IL1a_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/t2d_hm3/Pericytes_Adipose_basal_vs_IL1a_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/total_cholesterol/Pericytes_Adipose_basal_vs_IL1a_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/triacylglycerol/Pericytes_Adipose_basal_vs_IL1a_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/1707/Pericytes_Adipose_basal_vs_IL6TNFa_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/T1D/Pericytes_Adipose_basal_vs_IL6TNFa_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/aging/Pericytes_Adipose_basal_vs_IL6TNFa_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_vs_IL6TNFa_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/glaucoma/Pericytes_Adipose_basal_vs_IL6TNFa_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/ischemic_stroke/Pericytes_Adipose_basal_vs_IL6TNFa_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/t2d_2024/Pericytes_Adipose_basal_vs_IL6TNFa_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/t2d_hm3/Pericytes_Adipose_basal_vs_IL6TNFa_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/total_cholesterol/Pericytes_Adipose_basal_vs_IL6TNFa_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/triacylglycerol/Pericytes_Adipose_basal_vs_IL6TNFa_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/1707/Pericytes_Adipose_basal_vs_low_serum_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/T1D/Pericytes_Adipose_basal_vs_low_serum_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/aging/Pericytes_Adipose_basal_vs_low_serum_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_vs_low_serum_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/glaucoma/Pericytes_Adipose_basal_vs_low_serum_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/ischemic_stroke/Pericytes_Adipose_basal_vs_low_serum_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/t2d_2024/Pericytes_Adipose_basal_vs_low_serum_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/t2d_hm3/Pericytes_Adipose_basal_vs_low_serum_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/total_cholesterol/Pericytes_Adipose_basal_vs_low_serum_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/triacylglycerol/Pericytes_Adipose_basal_vs_low_serum_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/1707/Pericytes_Adipose_basal_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/T1D/Pericytes_Adipose_basal_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/aging/Pericytes_Adipose_basal_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/glaucoma/Pericytes_Adipose_basal_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/ischemic_stroke/Pericytes_Adipose_basal_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/t2d_2024/Pericytes_Adipose_basal_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/t2d_hm3/Pericytes_Adipose_basal_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/total_cholesterol/Pericytes_Adipose_basal_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/triacylglycerol/Pericytes_Adipose_basal_vs_oxLDL_triacylglycerol.results',  
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/1707/Pericytes_Adipose_low_serum_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/T1D/Pericytes_Adipose_low_serum_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/aging/Pericytes_Adipose_low_serum_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_low_serum_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/glaucoma/Pericytes_Adipose_low_serum_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/ischemic_stroke/Pericytes_Adipose_low_serum_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/t2d_2024/Pericytes_Adipose_low_serum_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/t2d_hm3/Pericytes_Adipose_low_serum_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/total_cholesterol/Pericytes_Adipose_low_serum_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/triacylglycerol/Pericytes_Adipose_low_serum_vs_oxLDL_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/1707/Pericytes_Adipose_basal_tstat_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/T1D/Pericytes_Adipose_basal_tstat_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/aging/Pericytes_Adipose_basal_tstat_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/CAD_GWAS_BBJ_meta_hm3/Pericytes_Adipose_basal_tstat_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/glaucoma/Pericytes_Adipose_basal_tstat_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/ischemic_stroke/Pericytes_Adipose_basal_tstat_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/t2d_2024/Pericytes_Adipose_basal_tstat_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/t2d_hm3/Pericytes_Adipose_basal_tstat_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/total_cholesterol/Pericytes_Adipose_basal_tstat_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/triacylglycerol/Pericytes_Adipose_basal_tstat_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/6142/Pericytes_Adipose_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/Height/Pericytes_Adipose_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/PA/basal_tstat/03_h2/SCZ/Pericytes_Adipose_basal_tstat_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/6142/Pericytes_Adipose_basal_vs_hg_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/6142/Pericytes_Adipose_basal_vs_IL1a_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/6142/Pericytes_Adipose_basal_vs_IL6TNFa_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/6142/Pericytes_Adipose_basal_vs_low_serum_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/6142/Pericytes_Adipose_basal_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/6142/Pericytes_Adipose_low_serum_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/Height/Pericytes_Adipose_basal_vs_hg_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/Height/Pericytes_Adipose_basal_vs_IL1a_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/Height/Pericytes_Adipose_basal_vs_IL6TNFa_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/Height/Pericytes_Adipose_basal_vs_oxLDL_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/Height/Pericytes_Adipose_basal_vs_low_serum_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_low_serum/03_h2/SCZ/Pericytes_Adipose_basal_vs_low_serum_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_hg/03_h2/SCZ/Pericytes_Adipose_basal_vs_hg_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL1a/03_h2/SCZ/Pericytes_Adipose_basal_vs_IL1a_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_IL6TNFa/03_h2/SCZ/Pericytes_Adipose_basal_vs_IL6TNFa_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/basal_vs_oxLDL/03_h2/SCZ/Pericytes_Adipose_basal_vs_oxLDL_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/Height/Pericytes_Adipose_low_serum_vs_oxLDL_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/PA/low_serum_vs_oxLDL/03_h2/SCZ/Pericytes_Adipose_low_serum_vs_oxLDL_SCZ.results'
)

vma_paths <- c(
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/1707/VSMC_pulmonary_artery_basal_vs_hg_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/T1D/VSMC_pulmonary_artery_basal_vs_hg_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/aging/VSMC_pulmonary_artery_basal_vs_hg_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_vs_hg_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/glaucoma/VSMC_pulmonary_artery_basal_vs_hg_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_vs_hg_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_vs_hg_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_vs_hg_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_vs_hg_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_vs_hg_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/1707/VSMC_pulmonary_artery_basal_vs_IL1a_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/T1D/VSMC_pulmonary_artery_basal_vs_IL1a_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/aging/VSMC_pulmonary_artery_basal_vs_IL1a_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_vs_IL1a_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/glaucoma/VSMC_pulmonary_artery_basal_vs_IL1a_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_vs_IL1a_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_vs_IL1a_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_vs_IL1a_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_vs_IL1a_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_vs_IL1a_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/1707/VSMC_pulmonary_artery_basal_vs_IL6TNFa_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/T1D/VSMC_pulmonary_artery_basal_vs_IL6TNFa_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/aging/VSMC_pulmonary_artery_basal_vs_IL6TNFa_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_vs_IL6TNFa_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/glaucoma/VSMC_pulmonary_artery_basal_vs_IL6TNFa_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_vs_IL6TNFa_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_vs_IL6TNFa_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_vs_IL6TNFa_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_vs_IL6TNFa_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_vs_IL6TNFa_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/1707/VSMC_pulmonary_artery_basal_vs_low_serum_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/T1D/VSMC_pulmonary_artery_basal_vs_low_serum_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/aging/VSMC_pulmonary_artery_basal_vs_low_serum_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_vs_low_serum_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/glaucoma/VSMC_pulmonary_artery_basal_vs_low_serum_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_vs_low_serum_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_vs_low_serum_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_vs_low_serum_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_vs_low_serum_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_vs_low_serum_triacylglycerol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/1707/VSMC_pulmonary_artery_basal_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/T1D/VSMC_pulmonary_artery_basal_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/aging/VSMC_pulmonary_artery_basal_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/glaucoma/VSMC_pulmonary_artery_basal_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_vs_oxLDL_triacylglycerol.results',  
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/1707/VSMC_pulmonary_artery_low_serum_vs_oxLDL_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/T1D/VSMC_pulmonary_artery_low_serum_vs_oxLDL_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/aging/VSMC_pulmonary_artery_low_serum_vs_oxLDL_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_low_serum_vs_oxLDL_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/glaucoma/VSMC_pulmonary_artery_low_serum_vs_oxLDL_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/ischemic_stroke/VSMC_pulmonary_artery_low_serum_vs_oxLDL_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/t2d_2024/VSMC_pulmonary_artery_low_serum_vs_oxLDL_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/t2d_hm3/VSMC_pulmonary_artery_low_serum_vs_oxLDL_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/total_cholesterol/VSMC_pulmonary_artery_low_serum_vs_oxLDL_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/triacylglycerol/VSMC_pulmonary_artery_low_serum_vs_oxLDL_triacylglycerol.results',  
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/1707/VSMC_pulmonary_artery_basal_tstat_1707.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/T1D/VSMC_pulmonary_artery_basal_tstat_T1D.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/aging/VSMC_pulmonary_artery_basal_tstat_aging.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/CAD_GWAS_BBJ_meta_hm3/VSMC_pulmonary_artery_basal_tstat_CAD_GWAS_BBJ_meta_hm3.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/glaucoma/VSMC_pulmonary_artery_basal_tstat_glaucoma.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/ischemic_stroke/VSMC_pulmonary_artery_basal_tstat_ischemic_stroke.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/t2d_2024/VSMC_pulmonary_artery_basal_tstat_t2d_2024.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/t2d_hm3/VSMC_pulmonary_artery_basal_tstat_t2d_hm3.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/total_cholesterol/VSMC_pulmonary_artery_basal_tstat_total_cholesterol.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/triacylglycerol/VSMC_pulmonary_artery_basal_tstat_triacylglycerol.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/6142/VSMC_pulmonary_artery_basal_tstat_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/Height/VSMC_pulmonary_artery_basal_tstat_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_tstats/VMA/basal_tstat/03_h2/SCZ/VSMC_pulmonary_artery_basal_tstat_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/6142/VSMC_pulmonary_artery_basal_vs_hg_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/6142/VSMC_pulmonary_artery_basal_vs_IL1a_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/6142/VSMC_pulmonary_artery_basal_vs_IL6TNFa_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/6142/VSMC_pulmonary_artery_basal_vs_low_serum_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/6142/VSMC_pulmonary_artery_basal_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/6142/VSMC_pulmonary_artery_low_serum_vs_oxLDL_6142.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/Height/VSMC_pulmonary_artery_basal_vs_hg_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/Height/VSMC_pulmonary_artery_basal_vs_IL1a_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/Height/VSMC_pulmonary_artery_basal_vs_IL6TNFa_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/Height/VSMC_pulmonary_artery_basal_vs_oxLDL_Height.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/Height/VSMC_pulmonary_artery_basal_vs_low_serum_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_low_serum/03_h2/SCZ/VSMC_pulmonary_artery_basal_vs_low_serum_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_hg/03_h2/SCZ/VSMC_pulmonary_artery_basal_vs_hg_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL1a/03_h2/SCZ/VSMC_pulmonary_artery_basal_vs_IL1a_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_IL6TNFa/03_h2/SCZ/VSMC_pulmonary_artery_basal_vs_IL6TNFa_SCZ.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/basal_vs_oxLDL/03_h2/SCZ/VSMC_pulmonary_artery_basal_vs_oxLDL_SCZ.results',
  # '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/Height/VSMC_pulmonary_artery_low_serum_vs_oxLDL_Height.results',
  '/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/VMA/low_serum_vs_oxLDL/03_h2/SCZ/VSMC_pulmonary_artery_low_serum_vs_oxLDL_SCZ.results'
)




# Define cell_type_paths 
cell_type_paths <- list(
  Cardiac_Fibroblasts_Atrial = cfa_paths,
  Cardiac_Fibroblasts_Ventricular = cfv_paths,
  Endothelial_Coronary_Artery = ec_paths,
  Pericytes_Adipose = pa_paths,
  VSMC_Pulmonary_Artery = vma_paths
)

# Function to load and combine results from file paths
load_and_combine_results <- function(file_paths, cell_type) {
  combined_results <- lapply(file_paths, function(path) {
    if (file.exists(path)) {
      results <- fread(path)
      condition <- str_extract(path, 
                               paste(c("basal_vs_hg", "basal_vs_IL1a", "basal_vs_IL6TNFa", 
                                       "basal_vs_low_serum", "basal_vs_oxLDL", "low_serum_vs_oxLDL", "basal"), 
                                     collapse = "|"))
      gwas <- str_extract(path, 
                          paste(c("SCZ", "1707", "CAD_GWAS_BBJ_meta_hm3",
                                  "ischemic_stroke", "t2d_2024", "total_cholesterol", "triacylglycerol"), 
                                collapse = "|"))
      results <- results %>%
        filter(Category == "L2_1") %>%
        mutate(
          GWAS = gwas,
          Condition = condition,
          Cell_Type = cell_type
        )
      return(results)
    } else {
      warning(sprintf("File not found: %s", path))
      return(NULL)
    }
  })
  
  combined_results <- combined_results[!sapply(combined_results, is.null)] %>% bind_rows()
  if (nrow(combined_results) == 0) {
    warning(sprintf("No valid data for cell type: %s", cell_type))
    return(data.frame())
  }
  combined_results <- combined_results %>%
    mutate(
      Condition = factor(Condition, levels = c("basal_vs_hg", "basal_vs_IL1a", "basal_vs_IL6TNFa", 
                                               "basal_vs_low_serum", "basal_vs_oxLDL", "low_serum_vs_oxLDL","basal")),
      GWAS = factor(GWAS, levels = c("CAD_GWAS_BBJ_meta_hm3", 
                                     "ischemic_stroke", "t2d_2024", "total_cholesterol", 
                                     "triacylglycerol", "1707", "SCZ")),
      Cell_Type = factor(Cell_Type)
    )
  return(combined_results)
}

# Load and combine results for all cell types
all_results <- lapply(names(cell_type_paths), function(cell_type) {
  file_paths <- cell_type_paths[[cell_type]]
  load_and_combine_results(file_paths, cell_type)
}) %>% bind_rows()

# Define custom GWAS trait labels
new_trait_labels <- c(
  # "T1D" = "T1D",
  # "aging" = "Aging",
  "CAD_GWAS_BBJ_meta_hm3" = "CAD",
  # "glaucoma" = "Glaucoma",
  "ischemic_stroke" = "Ischemic Stroke",
  "t2d_2024" = "T2D",
  "total_cholesterol" = "Total Cholesterol", 
  "triacylglycerol" = "Triglycerides",
  "1707" = "Left Handedness",
  # "6142" = "Employment Status",
  # "Height" = "Height",
  "SCZ" = "Schizophrenia"
)

# Function to map GWAS labels
map_gwas_labels <- function(gwas) {
  if (gwas %in% names(new_trait_labels)) {
    return(new_trait_labels[gwas])
  } else {
    return(gsub("_", " ", gwas))
  }
}

all_results <- all_results %>%
  mutate(GWAS_Label = sapply(as.character(GWAS), map_gwas_labels))

all_results <- all_results %>%
  mutate(GWAS_Label = factor(GWAS_Label,
                             levels = c("CAD", "Ischemic Stroke", "Total Cholesterol", "Triglycerides", "T2D", "Left Handedness", "Schizophrenia")))

all_results <- all_results %>%
  mutate(Condition = factor(Condition,
                            levels = c("low_serum_vs_oxLDL",
                                       "basal_vs_oxLDL", 
                                       "basal_vs_low_serum", 
                                       "basal_vs_hg",
                                       "basal_vs_IL6TNFa", 
                                       "basal_vs_IL1a",
                                       "basal")))

print(unique(all_results$GWAS))
print(unique(all_results$GWAS_Label))

plot_combined_dot_plot <- function(combined_results, significance_threshold = 3.35, n_categories = 112) {
  
  if (n_categories <= 0) {
    stop("Number of categories (n_categories) must be greater than zero.")
  }
  
  filtered_results <- combined_results %>%
    filter(!is.na(Enrichment), !is.na(Enrichment_p)) %>%
    mutate(
      # Apply Bonferroni correction: adjusted_P = raw p-value / number of tests
      adjusted_P = ifelse(Enrichment_p > 0, Enrichment_p / n_categories, .Machine$double.xmin),  # Smallest nonzero value if p = 0
      # Compute -log10 of adjusted p-value, preventing log10(0) issues
      log10_adjusted_P = ifelse(adjusted_P > 0, -log10(adjusted_P), NA_real_),
      # Define significance based on adjusted threshold
      Significant = log10_adjusted_P > significance_threshold,
      color_intensity = log10_adjusted_P,
      dot_size = Enrichment
    )
  
  print(paste("Number of rows after filtering:", nrow(filtered_results)))
  print(summary(filtered_results$Enrichment_p))  # Check p-value distribution
  
  if (nrow(filtered_results) == 0) {
    warning("No valid data after filtering. Check your input dataset.")
    return(NULL)
  }
  
  ggplot(filtered_results, aes(x = GWAS_Label, y = Condition)) +
    geom_point(aes(size = dot_size, fill = color_intensity),
               shape = 21,
               color = ifelse(filtered_results$Significant, "black", "transparent"),
               stroke = 0.5) + 
    scale_fill_gradient(low = "white", high = "red", name = "-log10(Adj P)",
                        guide = guide_colorbar(barwidth = 1, barheight = 10)) +
    scale_size_continuous(name = "Enrichment\nCoefficient", range = c(1, 8)) +
    labs(x = 'GWAS Traits', y = 'Condition') +
    facet_wrap(~ Cell_Type, scales = "free_y", 
               labeller = labeller(Cell_Type = function(x) gsub("_", " ", x))) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey80"),
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', colour = "grey90"),
      strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 1),
      strip.text = element_text(size = 14, face = "bold", colour = "black")
    )
}

# Generate the plot with improved Bonferroni correction
final_plot <- plot_combined_dot_plot(all_results, significance_threshold = 3.35, n_categories = 112)
final_plot
  
# Save the plot 
ggsave("/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA_Heritability_Enrichment_BC_w_tstat_basal.jpg",
         plot = final_plot, device = "jpeg", height = 5.5, width = 6.5, units = "in", dpi = 900)
ggsave("/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA_Heritability_Enrichment_BC_w_tstat_basal.pdf",
       plot = final_plot, device = "pdf", height = 5.5, width = 6.5, units = "in")
ggsave("/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/CFA_Heritability_Enrichment_BC_w_tstat_basal.svg",
       plot = final_plot, device = "svg", height = 5.5, width = 6.5, units = "in")




# ---- Generate and Save Supplementary Table ----

library(openxlsx)
write_supplementary_table <- function(results_df, output_path_csv, output_path_excel, n_categories = 112) {
  supp_table <- results_df %>%
    filter(!is.na(Enrichment), !is.na(Enrichment_p)) %>%
    mutate(`Adjusted P-Value` = pmin(Enrichment_p / n_categories, 1.0)) %>%
    select(
      `Cell Type` = Cell_Type,
      `Condition` = Condition,
      `GWAS Trait` = GWAS_Label,
      `Enrichment Coefficient` = Enrichment,
      `Nominal P-Value` = Enrichment_p,
      `Bonferroni Adjusted P-Value` = `Adjusted P-Value`
    ) %>%
    arrange(`Cell Type`, `Condition`, `GWAS Trait`)
  
  if (!is.null(output_path_csv)) {
    fwrite(supp_table, output_path_csv)
    message(paste("Supplementary table saved to CSV:", output_path_csv))
  }
  
  if (!is.null(output_path_excel)) {
    wb <- createWorkbook()
    addWorksheet(wb, "Heritability_Enrichment")
    writeData(wb, "Heritability_Enrichment", supp_table, withFilter = TRUE)
    setColWidths(wb, "Heritability_Enrichment", cols = 1:ncol(supp_table), widths = "auto")
    saveWorkbook(wb, output_path_excel, overwrite = TRUE)
    message(paste("Supplementary table saved to Excel:", output_path_excel))
  }
}

supp_table_csv_path <- "/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/Heritability_Enrichment_Supp_Table.csv"
supp_table_excel_path <- "/Volumes/broad_mcl/members_dir/mmurali/projects/9p21/terra_results_deg_hg38/Heritability_Enrichment_Supp_Table.xlsx"

write_supplementary_table(
  results_df = all_results,
  output_path_csv = supp_table_csv_path,
  output_path_excel = supp_table_excel_path,
  n_categories = 112
)
