rm(list = ls())
library(dplyr)
library(ggplot2)
library(tibble)

get_feature_names <- function(all) {
  all %>%
    mutate(feature_measure = case_when(
      grepl("Texture",            features) ~ "Texture",
      grepl("Intensity",          features) ~ "Intensity",
      grepl("Granularity",        features) ~ "Granularity",
      grepl("Correlation",        features) ~ "Correlation",
      grepl("RadialDistribution", features) ~ "RadialDistribution",
      TRUE                                  ~ "Location/Shape/Count"
    ))
}

add_channel_names <- function(all) {
  all %>%
    mutate(feature_channel = case_when(
      grepl("SmallBODIPY|LargeBODIPY", features)                               ~ "BODIPYRNA",
      grepl("Mito",      features) & !grepl("AGP|DNA|BODIPYRNA", features)     ~ "Mito",
      grepl("BODIPYRNA", features) & !grepl("Mito|AGP|DNA", features)          ~ "BODIPYRNA",
      grepl("AGP",       features) & !grepl("Mito|DNA|BODIPYRNA", features)    ~ "AGP",
      grepl("DNA",       features) & !grepl("Mito|AGP|BODIPYRNA", features)    ~ "DNA",
      rowSums(cbind(
        grepl("Mito",      features),
        grepl("AGP",       features),
        grepl("DNA",       features),
        grepl("BODIPYRNA", features)
      )) > 1                                                                   ~ "Combination",
      TRUE                                                                     ~ "Compartment"
    ))
}

channel_color <- c(
  AGP         = "#F2CF40",
  DNA         = "#65A7DB",
  Lipid       = "#4FAC2C",
  BODIPYRNA   = "#4FAC2C",
  Mito        = "#F56364",
  Compartment = "#A2A79F",
  Combination = "#818C7D"
)

get_channel_colors <- function(all) {
  all %>% mutate(color = channel_color[feature_channel])
}


#-------------Effect of KD
cfb <- read.csv("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/cfb_KD.csv") %>%
  filter(qvalue <= 0.2) %>%
  rename(features = featureName) %>%
  get_feature_names() %>%
  add_channel_names() %>%
  get_channel_colors()

smc <- read.csv("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/smc_KD.csv") %>%
  filter(qvalue <= 0.2) %>%
  rename(features = featureName) %>%
  get_feature_names() %>%
  add_channel_names() %>%
  get_channel_colors()

combined <- bind_rows(
  mutate(cfb, cell = "CFB\n qvalue ≤ 0.2"),
  mutate(smc, cell = "SMC\n qvalue ≤ 0.2")
)
write.csv(
  combined,
  "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/KD_featureCategorization.csv",
  row.names = FALSE
)

combined <- combined %>%
  mutate(
    stimulation = case_when(
      stimulation == "basal" ~ "Non-targeting control\nvs\nsiMTAP\n(baseline)",
      stimulation == "tgfb"  ~ "Non-targeting control\nvs\nsiMTAP\n(TGFb)",
      TRUE ~ stimulation
    )
  )

combined$stimulation <- factor(combined$stimulation, levels = unique(combined$stimulation))

y_max <- combined %>%
  count(cell, stimulation, feature_channel, name = "feature_count") %>%
  group_by(cell, stimulation) %>%
  summarise(total_count = sum(feature_count), .groups = "drop") %>%
  pull(total_count) %>%
  max() + 20

cfb_count <- combined %>%
  filter(cell == "CFB\n qvalue ≤ 0.2") %>%
  count(stimulation, feature_channel, name = "feature_count")

smc_count <- combined %>%
  filter(cell == "SMC\n qvalue ≤ 0.2") %>%
  count(stimulation, feature_channel, name = "feature_count")

p_cfb <- cfb_count %>%
  ggplot(aes(x = stimulation, y = feature_count, fill = feature_channel)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = filter(cfb_count, feature_count > 0),
    aes(label = feature_count),
    position = position_stack(vjust = 0.5),
    size = 5,
    fontface = "bold"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = channel_color) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
  labs(
    title = "CFB",
    x = "Stimulation",
    y = "Feature Count",
    fill = "Channel"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 18, face = "bold"),
    legend.title  = element_text(size = 16, face = "bold"),
    legend.text   = element_text(size = 14, face = "bold"),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

p_smc <- smc_count %>%
  ggplot(aes(x = stimulation, y = feature_count, fill = feature_channel)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = filter(smc_count, feature_count > 0),
    aes(label = feature_count),
    position = position_stack(vjust = 0.5),
    size = 6,
    fontface = "bold"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = channel_color) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
  labs(
    title = "SMC",
    x = "Stimulation",
    y = "Feature Count",
    fill = "Channel"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text     = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 18, face = "bold"),
    legend.title  = element_text(size = 16, face = "bold"),
    legend.text   = element_text(size = 14, face = "bold"),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

ggsave(
  "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/KD_cfb.png",
  plot = p_cfb, width = 10, height = 6
)
ggsave(
  "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/KD_cfb.svg",
  plot = p_cfb, width = 10, height = 6
)
ggsave(
  "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/KD_smc.png",
  plot = p_smc, width = 10, height = 6
)
ggsave(
  "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/KD_smc.svg",
  plot = p_smc, width = 10, height = 6
)

dev.off()
print("done")



#-------------Effect of Stim
# 
# cfbs <- read.csv("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/cfb_stim.csv") %>%
#   filter(qvalue <= 0.2) %>%
#   rename(features = featureName) %>%
#   get_feature_names() %>%
#   add_channel_names() %>%
#   get_channel_colors()
# 
# smcs <- read.csv("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/smc_stim.csv") %>%
#   filter(qvalue <= 0.2) %>%
#   rename(features = featureName) %>%
#   get_feature_names() %>%
#   add_channel_names() %>%
#   get_channel_colors()
# 
# cfbs <- cfbs %>%
#   mutate(basal = ifelse(basal == "nontargeting control", "nontargeting\ncontrol", basal),
#          cell = "CFB\n qvalue ≤ 0.2")
# 
# smcs <- smcs %>%
#   mutate(basal = ifelse(basal == "nontargeting control", "nontargeting\ncontrol", basal),
#          cell = "SMC\n qvalue ≤ 0.2")
# 
# combined_stim <- bind_rows(cfbs, smcs)
# 
# head(combined_stim)
# write.csv(combined_stim, "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/stim_featureCategorization.csv", row.names = FALSE)
# 
# y_max_stim <- combined_stim %>%
#   count(cell, basal, feature_channel, name = "feature_count") %>%
#   group_by(cell, basal) %>%
#   summarise(total_count = sum(feature_count)) %>%
#   pull(total_count) %>%
#   max() + 20
# 
# p<- combined_stim %>%
#   count(cell, basal, feature_channel, name = "feature_count") %>%
#   ggplot(aes(x = basal, y = feature_count, fill = feature_channel)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ cell) +
#   scale_fill_manual(values = channel_color) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, y_max_stim)) +
#   labs(
#     x = "Knockdown",
#     y = "Feature Count",
#     fill = "Channel"
#   ) +
#   theme_minimal(base_size = 16) +
#   theme(
#     axis.text = element_text(size = 16, face = "bold"),
#     axis.title = element_text(size = 18, face = "bold"),
#     legend.title = element_text(size = 16, face = "bold"),
#     legend.text = element_text(size = 14, face = "bold"),
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
#     panel.grid.major = element_line(color = "grey80"),
#     panel.grid.minor = element_line(color = "grey90"),
#     strip.text = element_text(size = 16, face = "bold")
#   )
# 
# ggsave("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/stim.png", plot = p, width = 10, height = 6)
# ggsave("/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/stim.svg", plot = p, width = 10, height = 6)
# 
# 
# 
# print("done")
































# uf <- unique(df$featureName)
# write.csv(data.frame(featureName = uf), 
#           "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/wilcoxAndStack/uf.csv", 
#           row.names = FALSE)

