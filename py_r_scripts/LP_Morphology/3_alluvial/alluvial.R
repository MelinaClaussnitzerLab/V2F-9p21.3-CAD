library(googleCloudStorageR)
library(readr)
library(googledrive)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(moments)
library(RNOmni)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(coin)
Sys.setenv(VROOM_CONNECTION_SIZE = 10485760)
library(ggalluvial)   
library(qvalue)
library(alluvial)
suppressPackageStartupMessages( require(tidyverse) )
suppressPackageStartupMessages( require(easyalluvial) )
library(umap)
library(viridis)
library(colorspace) 
library(grDevices)
library(RColorBrewer)
library(randomcoloR)
library(cowplot)
library(grid)
options(readr.show_col_types = FALSE)
options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)
library(knitr)
library(gtools)
library(pheatmap)
library(scales)
library(knitr)
library(kableExtra)
library(grDevices)
library(sp)
library(ggdendro)
library(ggfortify)
library(ggrepel)
# library(WMDB)
# library(NMF)
library(plotly)
library(pheatmap)


raw_df <- read_csv("/home/jupyter/9p21_071825/raw_df.csv")
head(raw_df)



metaindex <- grep("Metadata|Condition|Batch|Plate_Map|color|Image", colnames(raw_df))
features <- colnames(raw_df)[-metaindex]

all <- data.frame(features)
colnames(all) <- "features"
all$features <- gsub("LP_", "", all$features)


unique_parts <- unique(sapply(all$features, function(x) strsplit(x, "_")[[1]][1]))



all <- all %>%
  mutate(Compartment = case_when(
    grepl("Cells_", features) ~ "Cells",
    grepl("Cytoplasm_", features) ~ "Cytoplasm",
    grepl("Nuclei", features) ~ "Nuclei"
  ))


all <- all %>%
  mutate(Channel = case_when(
    grepl("SmallBODIPY|LargeBODIPY", features) ~ "BODIPYRNA",
    
    grepl("Mito", features) &
      !grepl("AGP|BODIPYRNA|DNA|SmallBODIPY|LargeBODIPY", features) ~ "Mito",
    
    grepl("BODIPYRNA", features) &
      !grepl("AGP|Mito|DNA", features) ~ "BODIPYRNA",
    
    grepl("AGP", features) &
      !grepl("Mito|BODIPYRNA|DNA|SmallBODIPY|LargeBODIPY", features) ~ "AGP",
    
    grepl("DNA", features) &
      !grepl("AGP|Mito|BODIPYRNA|SmallBODIPY|LargeBODIPY", features) ~ "DNA",
    
    (rowSums(cbind(
      grepl("Mito", features),
      grepl("AGP", features),
      grepl("DNA", features),
      grepl("BODIPYRNA", features)
    )) > 1) ~ "Combination",
    
    TRUE ~ "Compartment"
  ))

all <- all %>% 
  mutate(Measure = case_when(
    grepl("Intensity", features) &
      !grepl("Texture|Granularity|Correlation|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", features) 
    ~ "Intensity",
    
    grepl("Texture", features) &
      !grepl("Intensity|Granularity|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", features) 
    ~ "Texture",
    
    grepl("Granularity", features) &
      !grepl("Intensity|Texture|Correlation|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", features) 
    ~ "Granularity",
    
    grepl("Correlation", features) &
      !grepl("Texture|Intensity|Granularity|AreaShape|RadialDistribution|Location|Neighbors|Count|Number", features) 
    ~ "Correlation",
    
    grepl("RadialDistribution", features) &
      !grepl("Intensity|Texture|Granularity|Correlation|AreaShape|Location|Neighbors|Count|Number", features) 
    ~ "RadialDistribution",
    
    TRUE 
    ~ "Location/Shape/Size/\nNeighbor/Count"
  ))


col_vector <- c(
  "#A2A89F", "#55828B", "#87BBA2", "#f2cf41", "#50ac2c",
  "#828C7D", "#A2A89F", "#65a6db", "#f56464", "#828C7D",
  "#C9E4CA", "#A1D685"
)
col_vector_compartments <- c("#836953", "#306178", "#8965b3")

p <- alluvial_wide(
  select(all, Compartment, Channel, Measure),
  fill_by          = "first_variable",
  order_levels     = as.character(col_vector),
  col_vector_value = col_vector,
  col_vector_flow  = col_vector_compartments,
  stratum_width    = 0.25
)

p$layers[[3]] <- NULL

gb       <- ggplot_build(p)
strataDF <- gb$data[[2]]

p_final <- p +
  geom_text(
    data        = strataDF,
    aes(x = x, y = y, label = stratum),
    size        = 5,
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(
    panel.background = element_blank(),
    text             = element_text(face = "bold"),
    axis.text.x      = element_blank(),
    axis.title.x     = element_blank(),
    axis.title.y = element_text(size = 70, face = "bold"),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 50, face = "bold"),
    plot.caption     = element_blank(),
    plot.margin      = unit(c(1, 3, 1, 1), "inches"),   
    panel.grid       = element_blank()
  )

print(p_final)






































