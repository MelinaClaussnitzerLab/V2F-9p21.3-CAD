library(readr)
library(dplyr)
library(RNOmni)
library(ggplot2)
library(tidyr)
df <- read.csv('/Users/cho/Documents/cellprofiler/downstream/9p21/072225/data/data_.csv')

df <- df %>% filter(Metadata_cell_type == 'primary cfb')



print(any(is.na(df)))

meta_index <- grep("Metadata|Condition|Batch|Plate_Map|Image", colnames(df))
meta <- dplyr::select(df, all_of(meta_index))
colnames(meta)


#-------drop 0 std dev
zero_sd_cols <- names(df)[
  sapply(df, function(col) is.numeric(col) && sd(col, na.rm = TRUE) == 0) &
    !grepl("Metadata", names(df))
]

cat("Removing", length(zero_sd_cols), "columns with zero standard deviation (excluding Metadata):\n")

if (length(zero_sd_cols) > 0) {
  for (col in zero_sd_cols) {
    vals <- unique(df[[col]])
    cat(" •", col, "->", paste(vals, collapse = ", "), "\n")
  }
} else {
  cat(" (none)\n")
}

df <- df %>% 
  select(-all_of(zero_sd_cols))



#--------normalize


lp <- dplyr::select(df, -all_of(meta_index))

normalized_data <- as.data.frame(apply(lp, 2, RankNorm))

norm_df <- cbind(meta, normalized_data)
head(norm_df)
print(any(is.na(norm_df)))

#----------save

write.csv(norm_df,"/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wellClustering/int_cfb.csv", row.names = FALSE)








#-------------hist
# orig_area <- df$Cells_AreaShape_Solidity
# int_area  <- norm_df$Cells_AreaShape_Solidity
# 
# plot_df <- data.frame(
#   Area    = c(orig_area, int_area),
#   Dataset = rep(c("Original", "INT"), each = length(orig_area))
# )
# 
# print(range(plot_df$Area))
# 
# library(ggplot2)
# ggplot(plot_df, aes(x = Area, fill = Dataset)) +
#   geom_density(alpha = 0.4, color = NA) +
#   scale_x_continuous(limits = range(plot_df$Area)) +
#   labs(
#     title = "Density of Cells_AreaShape_Area: Original vs. INT",
#     x     = "Cells_AreaShape_Area",
#     y     = "Density"
#   ) +
#   theme_minimal()



