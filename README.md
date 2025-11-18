# The <i>9p21.3</i> Coronary Artery Disease Risk Locus Modulates Vascular Cell-State Transitions via Enhancer-Driven Regulation of <i>MTAP</i>


## metadata
- Haplotypes:
  - [Haplotypes](Haplotypes.csv)
- FinnGen PheWAS:
  - https://r13.finngen.fi/variant/9:22098575-A-G


## Scripts
- ChromHMM scripts
  - Filter MintChIP peaks [filter_mintchip_macs_peaks.py](py_r_scripts/ChromHMM/000_filter_mintchip_macs_peaks.py)
  - Filter ATAC peaks [filter_atac_macs_peaks.py](py_r_scripts/ChromHMM/001_filter_atac_macs_peaks.py)
  - Jointly binerize [jointly_binarize_peaks.py](py_r_scripts/ChromHMM/002_jointly_binarize_peaks.py)
  - Jointly model [jointly_model.py](py_r_scripts/ChromHMM/003_jointly_model.py)

- Stratified LDSC
  - [readme](py_r_scripts/LDSC_heritability/README.md)
  - [analysis](py_r_scripts/LDSC_heritability)

- MAGMA analysis
  - [readme](py_r_scripts/MAGMA_analysis/README.md)
  - [analysis](py_r_scripts/MAGMA_analysis)
  

- LP_Morphology
  - [readme](py_r_scripts/LP_Morphology/README.md)
  - [analysis](py_r_scripts/LP_Morphology)

- Identifying novel enhancers
  - Extracting RoadMaps enhancer regions [100_extract_roadmap_enhancers.py](py_r_scripts/postprocessin_enhancers/100_extract_roadmap_enhancers.py)
  - Extracting VascCellWalls enhancer regions [101_parse_chromHMM_JointModel_enhancers.py](py_r_scripts/postprocessin_enhancers/101_parse_chromHMM_JointModel_enhancers.py)
  - Merge and intersect peaks [102_merge_and_intersect_bed_files.py](py_r_scripts/postprocessin_enhancers/102_merge_and_intersect_bed_files.py)
  - Extract novel enhancer regions [103_extract_novel_chromhmm_enhancers.py](py_r_scripts/postprocessin_enhancers/103_extract_novel_chromhmm_enhancers.py)

- ROSE 
  - Prepare ROSE input [00_make_enhancer_gff_and_scripts.py](py_r_scripts/ROSE/00_make_enhancer_gff_and_scripts.py)
  

- Single-nucleus RNA-seq from human hearts
  - Adjusting and plotting expression profiles [10_expressions.Rmd](py_r_scripts/snRNA-seq_human_hearts/10_expressions.Rmd)

- MAC-Seq
  - [QC](py_r_scripts/MAC_Seq/plate_qc_scripts)
  - [Sliding window](py_r_scripts/MAC_Seq/macseq_sliding_window_stats.R) 
  - [CFB analysis](py_r_scripts/MAC_Seq/cfb_final_analysis)
  - [SMC analysis](py_r_scripts/MAC_Seq/smc_final_analysis)

