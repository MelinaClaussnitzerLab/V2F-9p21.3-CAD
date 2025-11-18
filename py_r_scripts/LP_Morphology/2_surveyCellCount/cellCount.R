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

df <- read.csv("/Users/cho/Desktop/9p21Manuscript/cellCount/raw_df.csv")
head(df[1:10])

natural_sort <- function(x) {
  as.character(mixedsort(as.character(x)))
}

calculate_outliers <- function(df) {
  Q1 <- quantile(df$Metadata_Object_Count, 0.25)
  Q3 <- quantile(df$Metadata_Object_Count, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  sum(df$Metadata_Object_Count < lower_bound | df$Metadata_Object_Count > upper_bound)
}

outlier_counts <- df %>%
  group_by(Metadata_cell_type) %>%
  summarise(Outliers = calculate_outliers(cur_data()),
            Total = n())

outlier_counts <- outlier_counts %>%
  mutate(Metadata_cell_type_Outliers = factor(paste0(Metadata_cell_type, " (Outliers: ", Outliers, "/", Total, ")"),
                                              levels = natural_sort(paste0(Metadata_cell_type, " (Outliers: ", Outliers, "/", Total, ")"))))

plot_df <- df %>%
  left_join(outlier_counts, by = "Metadata_cell_type")

plot_df <- plot_df %>%
  group_by(Metadata_cell_type) %>%
  mutate(
    Q1 = quantile(Metadata_Object_Count, 0.25),
    Q3 = quantile(Metadata_Object_Count, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = Metadata_Object_Count < lower_bound | Metadata_Object_Count > upper_bound
  ) %>%
  ungroup() %>%
  mutate(row_num = row_number())

alternate_position <- function(row_num) {
  if (row_num %% 2 == 0) {
    return(-0.3) 
  } else {
    return(1.3) 
  }
}

plot_df <- plot_df %>%
  mutate(hjust_pos = sapply(row_num, alternate_position))

p <- ggplot(plot_df, aes(x = Metadata_cell_type_Outliers, y = Metadata_Object_Count)) +
  geom_violin(fill = "white", alpha = 0.3, color = "black") + 
  geom_boxplot(outlier.colour = NA, fill = "grey", width = 0.3) + 
  geom_jitter(data = subset(plot_df, !is_outlier), aes(x = Metadata_cell_type_Outliers, y = Metadata_Object_Count, fill = Metadata_stimulation),
              colour = "black", shape = 21, size = 1.5, width = 0.1) + 
  geom_point(data = subset(plot_df, is_outlier), aes(x = Metadata_cell_type_Outliers, y = Metadata_Object_Count),
             colour = "black", shape = 21, size = 3, fill = "red") +  
  geom_text(data = subset(plot_df, is_outlier), aes(label = paste(Metadata_treatment, Metadata_stimulation, sep = ":"), x = Metadata_cell_type_Outliers, y = Metadata_Object_Count, hjust = hjust_pos),
            vjust = 0, size = 3, colour = "black") +
  facet_wrap(~ Metadata_cell_type_Outliers, scales = "free_x", ncol = 3, strip.position = "top") +
  labs(title = "Box Plot of Cell Count by Cell Type",
       x = "Plate",
       y = "Count of Objects",
       fill = "Stimulation") +  # Adding the fill legend label
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.text.x = element_text(size = 10),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_line(colour = "grey", size = 0.5), 
    panel.grid.minor.y = element_line(colour = "grey", size = 0.25)
  ) +
  scale_fill_manual(values = c("basal" = "green", "tgfb" = "blue")) +
  scale_y_continuous(limits = range(df$Metadata_Object_Count))

print(p)




# dev.off()
# system(paste("gsutil cp", output_file, gcs_path))
