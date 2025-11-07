#

library(dplyr)
library(ggplot2)
library(reghelper)
library(patchwork)
library(ggrepel)
library(svglite)
library(tidyverse)
library(ggpubr)
library(ggnewscale)

####set savepath and read in data####

screen_lims = xlim(22064500,22163000)

loc_lims = xlim(22064500,22132000)

crop_lims = xlim(22095000,22120000)

input_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/1_platemaps_resources"
save_path_figs = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/locus_plots"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/datasets"
save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/sliding_window"
save_path_snps = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/enhancer_and_snp_plots"
save_path_tss = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/tss"

export_fig2 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/screen_fig"
export_fig4 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/supp_fig4"
export_fig2 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/fig2"

####enhancer dotplot####

enhancer_effectsize <- read.csv(paste0(save_path_snps,"/cfb_enhancer_mean_effects.csv"))

enhancer_effectsize <- enhancer_effectsize %>% filter(Gene != "MIR31HG")

enhancer_effectsize$Enhancer <- factor(enhancer_effectsize$Enhancer, levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12"),
                                       ordered = TRUE, labels=c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12"))

#enhancer_effectsize$Gene <- factor(enhancer_effectsize$Gene, levels = c("CDKN2B","MTAP","CDKN2A","CDKN2B-AS1","DMRTA1","MIR31HG"),
#                                       ordered = TRUE)

enhancer_effectsize$Gene <- factor(enhancer_effectsize$Gene, levels = c("DMRTA1","CDKN2B-AS1","CDKN2A","MTAP","CDKN2B"),
                                   ordered = TRUE)

#set legend for significance
color_key = c("not" = "white", "sig" = "black")
legend_order <- c("sig", "not")
size_order <- c("low","high")
size_key = c("low" = 8, "high" = 12)

enhancer_effectsize <- enhancer_effectsize %>% filter(Diff <= -0.1 | Diff >= 0.1)

enhancer_effectsize <- enhancer_effectsize %>% mutate(
  high_fold = case_when(abs(Diff) >= 0.2 ~ "high",
                        abs(Diff) < 0.2 ~ "low")
)

simpletheme <-   theme_bw() +
  theme(text = element_text(size = 20),
        title = element_text(size = 20),
        legend.text = element_text(size=10),
        axis.title.x=element_text(size = 25),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size = 25))

enhancer_effect_plot <- ggplot(data = enhancer_effectsize, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="blue", mid="white", high="red",
                         limits = c(-0.6, 0.6),
                         breaks = c(-0.6, -0.2, 0, 0.2, 0.6),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) + 
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "enhancer_dotplot.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_enhancer_dotplot.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)

#now remove all upregulated enhancers and repeat

enhancer_effectsize_down <- enhancer_effectsize %>% filter(Diff < -0.1)

enhancer_effect_plot <- ggplot(data = enhancer_effectsize_down, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("> -0.2", "\u2264 -0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="blue", high="white",
                         limits = c(-0.6,0),
                         breaks = c(-0.6,-0.3,0),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("> -0.2", "\u2264 -0.2")) +
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "enhancer_dotplot_downreg.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_enhancer_dotplot_downreg.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)

enhancer_effectsize_up <- enhancer_effectsize %>% filter(Diff >0.1)

enhancer_effect_plot <- ggplot(data = enhancer_effectsize_up, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="white", high="red",
                         limits = c(0,0.6),
                         breaks = c(0,0.3,0.6),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "enhancer_dotplot_upreg.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_enhancer_dotplot_upreg.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)

####PLOT INCLUDING MIR31HG####

enhancer_effectsize <- read.csv(paste0(save_path_snps,"/cfb_enhancer_mean_effects.csv"))
#enhancer_effectsize <- enhancer_effectsize %>% filter(Gene != "MIR31HG")
enhancer_effectsize$Enhancer <- factor(enhancer_effectsize$Enhancer, levels = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12"),
                                       ordered = TRUE, labels=c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12"))

enhancer_effectsize$Gene <- factor(enhancer_effectsize$Gene, levels = c("MIR31HG","DMRTA1","CDKN2B-AS1","CDKN2A","MTAP","CDKN2B"),
                                   ordered = TRUE)

enhancer_effectsize <- enhancer_effectsize %>% filter(Diff <= -0.1 | Diff >= 0.1)

enhancer_effectsize <- enhancer_effectsize %>% mutate(
  high_fold = case_when(abs(Diff) >= 0.2 ~ "high",
                        abs(Diff) < 0.2 ~ "low")
)

enhancer_effect_plot <- ggplot(data = enhancer_effectsize, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="blue", mid="white", high="red",
                         limits = c(-0.6, 0.6),
                         breaks = c(-0.6, -0.2, 0, 0.2, 0.6),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) + 
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "mir31hg_enhancer_dotplot.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_mir31hg_enhancer_dotplot.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)

#now remove all upregulated enhancers and repeat

enhancer_effectsize_down <- enhancer_effectsize %>% filter(Diff < -0.1)

enhancer_effect_plot <- ggplot(data = enhancer_effectsize_down, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="blue", high="white",
                         limits = c(-0.6,0),
                         breaks = c(-0.6,-0.3,0),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("> -0.2", "\u2264 -0.2")) +
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "mir31hg_enhancer_dotplot_downreg.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_mir31hg_enhancer_dotplot_downreg.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)


enhancer_effectsize_up <- enhancer_effectsize %>% filter(Diff >0.1)

enhancer_effect_plot <- ggplot(data = enhancer_effectsize_up, aes(x=Enhancer, y=Gene)) + 
  geom_point(aes(color=Diff, size= high_fold)) + #add: in aes, (size=FDR), foprmer: size = 12
  scale_radius() + scale_size_manual(breaks = size_order, values = size_key, labels = c("> -0.2", "\u2264 -0.2")) +
  #                                       limits = c(0,17),
  #                                       breaks = c(0,1,5,10,15)) +
  scale_y_discrete(drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradient2(low="white", high="red",
                         limits = c(0,0.6),
                         breaks = c(0,0.3,0.6),
                         #labels = c("-3", "-1", "0", "1", "3"),  
                         guide=guide_colorbar() )  +
  #set legend before changing color
  labs(colour="Gene expression\n(z-score)", title = "9p21 enhancer-gene interactions", size = "z-score magnitude") + #size="-log10 FDR", 
  xlab("Enhancer") + ylab("Gene") +
  #define new color scale for P value significance
  new_scale_color() +
  geom_point(aes(color = Signif_pcr, size = high_fold), shape = 1) +
  scale_color_manual(breaks = legend_order, values = color_key, guide = "none") + 
  scale_size_manual(breaks = size_order, values = size_key, labels = c("<0.2", "\u2265 0.2")) +
  #define rest of theme
  simpletheme
enhancer_effect_plot

ggsave(filename= file.path(save_path_snps, "mir31hg_enhancer_dotplot_upreg.png"), 
       plot=enhancer_effect_plot, device=png, width = 14, height = 4.2)

ggsave(filename= file.path(export_fig2, "cfb_mir31hg_enhancer_dotplot_upreg.svg"), 
       plot=enhancer_effect_plot, width = 14, height = 4.2)