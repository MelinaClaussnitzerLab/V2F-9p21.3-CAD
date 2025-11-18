f[9p21/LP_Morphology/0_cppipe/](#)

Description:
------------
This folder contains CellProfiler pipeline files related to image analysis and illumination correction for the 9p21 LP project. These pipelines are used to process high-content imaging data.

Contents:
---------
1. analysis.cppipe
   - This is the primary analysis pipeline. Segmentation and feature extraction happen here

2. illum.cppipe
   - This is the illumination correction pipeline.


[9p21/LP_Morphology/2_surveyCellCount/](#)

Description:
------------
This folder contains R scripts and raw data for visualizing and quantifying cell count distributions across conditions

Contents:
---------
1. cellCount.R  
   - R script that reads in the raw cell count data, detects outliers within each cell type group, and visualizes distributions using violin and box plots.  
   - The plot highlights outliers with red markers and labels them with treatment and stimulation metadata.  

2. raw_df.csv  
   - Raw input dataset used by `cellCount.R`.  
   - Includes metadata fields such as `Metadata_cell_type`, `Metadata_treatment`, `Metadata_stimulation`, and `Metadata_Object_Count`.  
   - This file provides the underlying measurements and annotations for cell count visualizations and analysis.

[9p21/LP_Morphology/3_alluvial/](#)

Description:
------------
This folder contains R scripts, raw data, and visual outputs related to the generation of alluvial plots summarizing feature compartment-channel-measure relationships

Contents:
---------
1. alluvial.R  
   - R script that processes feature names from `raw_df.csv` to classify them by cellular compartment (e.g., Nuclei, Cytoplasm), channel (e.g., DNA, AGP), and measurement type (e.g., Intensity, Texture).
   - Generates alluvial plots using `ggalluvial` to visualize how features distribute across compartments, channels, and measurement types.

2. Visualizations  
   - svg and pdf file outputs

[9p21/LP_Morphology/4_featureReduction/](#)

Description:
------------
This folder contains input, output, and clustering results from feature reduction

Contents:
---------
1. raw_df.csv  
   - Input dataset containing extracted features used for dimensionality reduction.

2. data_.csv  
   - Output data after performing feature reduction.

3. _clusters.csv  
   - Lists feature groupings based on clustering results, showing which features are grouped.

5. featureReduction.ipynb  
   - Jupyter notebook used to run feature reduction and clustering workflow.

[9p21/LP_Morphology/5_inverseNormalTransformation/](#)

Description:
------------
This folder contains inverse normal transformed data and results related to the normalization of morphological features for each cell type

Contents:
---------
1. InverseNormalTransformation.R  
   - R script used to apply inverse normal transformation to feature data

2. *_int.csv  
   - Inverse normal transformed feature data, saved separately by cell type:
     - `cfb_int.csv` – transformed data for CFB cells  
     - `smc_int.csv` – transformed data for SMC cells

3. data_.csv  
   - input of inverse normal transformation. output of feature reduction


[9p21/LP_Morphology/6_wilcoxTestAndStackPlot/](#)

Description:
------------
This folder contains scripts and data for statistical testing and visualization of siMTAP knockdown effects across cell types and stimulation conditions

Contents:
---------
1. *_KD.csv  
   - Results of testing the effect of siMTAP knockdown vs. control under each stimulation condition.  
   - These files serve as input to `stack.R` for downstream visualization.

2. wilcox.R  
   - Script that performs Wilcoxon rank-sum tests using inverse-normal transformed data.

3. stack.R  
   - Script for generating stacked plots based on knockdown test results.


[9p21/LP_Morphology/7_Morpheus_Clustering/](#)

Description:
------------
This folder contains input data and output visualizations from Broad Institute’s Morpheus software

Contents:
---------
1. *_INT_pearson_avg_input.csv  
   - Input matrices for clustering

2. *.pdf / *.svg  
   - Clustering heatmaps exported from Morpheus for both CFB and SMC datasets.

















