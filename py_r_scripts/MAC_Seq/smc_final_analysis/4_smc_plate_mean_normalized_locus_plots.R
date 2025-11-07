#Script for gaussian regression of MAC-Seq samples, split by cell type

library(dplyr)
library(ggplot2)
library(reghelper)
library(patchwork)
library(ggrepel)
library(svglite)
library(tidyverse)

####set savepath and read in data####

save_path = "output/glm/smc"
input_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/1_platemaps_resources"
save_path_figs = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/smc/zscore_normalized/locus_plots"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/smc/zscore_normalized/datasets"
save_path_batch = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/smc/zscore_normalized/batch_effect"
save_path_screen= "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/smc/zscore_normalized/whole_locus_vis"

screen_lims = xlim(22064500,22163000)

loc_lims = xlim(22064500,22132000)

crop_lims = xlim(22095000,22120000)

#read in smc data by replicate
smc_merged <- read.csv(paste0(save_path_files,"/plate_normalized_guide_means_output.csv"))
smc_merged <- smc_merged[,-1]

#convert metadata to factors for regression
smc_merged$sequence <- as.factor(smc_merged$sequence)
#smc_merged$cellType <- as.factor(smc_merged$cellType)
smc_merged$plate <- as.factor(smc_merged$plate)
smc_merged$rep_id <- as.factor(smc_merged$rep_id)

#define graphed position for negative controls (2 kinds), and assign these locations in the data
negguidepos = 22065000
locguidepos = 22070000

#plot histograms to check the normal distribution of these adjusted values
hist(smc_merged$CDKN2B.AS1, breaks=20, main="Normalized ANRIL expression per well", xlim=c(-15,15))
hist(smc_merged$CDKN2B, breaks=20, main="Normalized CDKN2B expression per well", xlim=c(-15,15))
hist(smc_merged$CDKN2A, breaks=20, main="Normalized CDKN2A expression per well" , xlim=c(-15,15))
hist(smc_merged$MTAP, breaks=20, main="Normalized MTAP expression per well", xlim=c(-15,15))
hist(smc_merged$DMRTA1, breaks=20, main="Normalized DMRTA1 expression per well", xlim=c(-15,15))
hist(smc_merged$MIR31HG, breaks=20, main="Normalized MIR31HG expression per well", xlim=c(-15,15))

#plot histograms to check the normal distribution of these adjusted values
hist(smc_merged$CDKN2B.AS1.mean, breaks=15, main="Normalized ANRIL expression per well", xlim=c(-15,15))
hist(smc_merged$CDKN2B.mean, breaks=15, main="Normalized CDKN2B expression per well", xlim=c(-15,15))
hist(smc_merged$CDKN2A.mean, breaks=15, main="Normalized CDKN2A expression per well" , xlim=c(-15,15))
hist(smc_merged$MTAP.mean, breaks=15, main="Normalized MTAP expression per well", xlim=c(-15,15))
hist(smc_merged$DMRTA1.mean, breaks=15, main="Normalized DMRTA1 expression per well", xlim=c(-15,15))
hist(smc_merged$MIR31HG.mean, breaks=20, main="Normalized MIR31HG expression per well", xlim=c(-15,15))


#save file of all readcounts and metadata for regression
#write.csv(smc_merged, paste0(save_path_files,"/plate_normalized_cpm_reads_all_input.csv"))

####make enhancer and SNP plots for use with other graphics####
#enhancer_key <- read.csv("raw_data/enhancer_key.csv")
enhancer_key <- read.csv(paste0(input_path,"/enhancer_key_controls.csv"))

stanford_key <- read.csv(paste0(input_path,"/enhancer_key_stanford_only.csv"))
enhancer_key_simple <- read.csv(paste0(input_path,"/enhancer_key.csv"))
enhancer_key_simple <- enhancer_key_simple %>% filter(Start < 22129200)

snp_key <- read.csv(paste0(input_path,"/SNP_key.csv"),fileEncoding="latin1")
snp_key$loc <- as.numeric(snp_key$loc)

#create simple enhancer graphic for combination with other plots
enhancer_plot <- ggplot(data=enhancer_key) +
  xlim(22064500,22132000) + #xlim(22064500,22163000) + ylim (-1,1) +
  geom_rect(data = enhancer_key, inherit.aes = F, aes(xmin = Start, xmax = End, ymin = -1,ymax = 0)) +
  geom_text_repel(data=enhancer_key, size=4, inherit.aes=F, aes(x=Mid, y=0.5, label=Enhancer),direction = "x") +
  theme_classic() + theme(
    text=element_text(size = 15),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size = 15),
    legend.position="none"
  ) + 
  labs(title = "9p21 enhancers", x="Chr.9")
enhancer_plot

#create simple enhancer graphic for combination with other plots
stanford_plot <- ggplot(data=stanford_key) +
  xlim(22064500,22132000) + #xlim(22064500,22163000) + ylim (-1,1) +
  #geom_rect(data = enhancer_key_simple, inherit.aes = F, aes(xmin = Start, xmax = End, ymin = -1,ymax = 0), alpha = 0.2) +
  geom_rect(data = stanford_key, inherit.aes = F, aes(xmin = Start, xmax = End, ymin = -1,ymax = 0)) +
  geom_text_repel(data=stanford_key, size=4, inherit.aes=F, aes(x=Mid, y=0.5, label=Enhancer),
                  direction = "x", box.padding = 0.2, force = 0.2) +
  theme_classic() + theme(
    text=element_text(size = 15),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size = 15),
    legend.position="none"
  ) + 
  labs(title = "9p21 Stanford enhancers", x="Chr.9")
stanford_plot

enhancer_plot_simple <- ggplot(data=enhancer_key_simple) +
  xlim(22064500,22132000) + #xlim(22064500,22163000) + ylim (-1,1) +
  #geom_rect(data = stanford_key, inherit.aes = F, aes(xmin = Start, xmax = End, ymin = -1,ymax = 0), alpha = 0.3) +
  geom_rect(data = enhancer_key_simple, inherit.aes = F, aes(xmin = Start, xmax = End, ymin = -1,ymax = 0)) +
  geom_text_repel(data=enhancer_key_simple, size=4, inherit.aes=F, aes(x=Mid, y=0.5, label=Enhancer),
                  direction = "x", box.padding = 0.2, force = 0.2) +
  theme_classic() + theme(
    text=element_text(size = 15),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size = 15),
    legend.position="none"
  ) + 
  labs(title = "9p21 Broad enhancers", x="Chr.9")
enhancer_plot_simple

temp <- enhancer_plot_simple/stanford_plot
ggsave(filename= file.path(save_path_figs, "broad_stanford_enhancer_compare.png"), plot=temp, device=png, width = 10, height = 4)
rm(temp)

#create simple SNP graphic for combination with other plots
snp_plot <- ggplot(data=snp_key) +
  geom_rect(data = snp_key, inherit.aes = F, aes(
    xmin = loc-10, xmax = loc+10, ymin = -0.05,ymax = 0.05, fill=is.na(constraint), color=is.na(constraint))) +
  xlim(22064500,22132000) + #xlim(22064500,22163000) + #ylim (-1,1) +
  #geom_text_repel(data=snp_key, size=4, inherit.aes=F, aes(x=GuidePosition, y=0.06, label=constraint)) +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "black"), name = "", labels = c("Constrained", "Not constrained")) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"), name = "", labels = c("Constrained", "Not constrained")) +
  #add enhancer rectangles with low alpha
  geom_rect(data = enhancer_key, inherit.aes = F, aes(
    xmin = Start, xmax = End, ymin = -Inf,ymax = Inf), alpha = 0.3, fill = "lightgray") +
  theme_classic() + theme(
    text=element_text(size = 15),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x=element_text(size = 15)
    #legend.position="none"
  ) + 
  labs(title = "9p21 SNPs", x="Chr.9")
snp_plot

####PLOT DATA ON LOCUS####

#subset(smc_merged, smc_merged$CDKN2A < (-0.99) & guide_ctrl == "Screen")
#  color="red") +
  
  #a=ggplot object
#b=coefs
#sd(smc_merged$CDKN2A)
nnfplot <- function(a,b,g) {
  y <- b %>% 
    summarize(mean=mean(b[[g]]),
              min=min(b[[g]]),
              max= max(b[[g]]),
              sd = sd(b[[g]]))
  c=as.numeric(y[1])
  d = as.numeric(y[2])
  e = as.numeric(y[3])
  f = as.numeric(y[4])
  h <- b %>% subset(b[[g]] < -1 & guide_ctrl == "Screen")
  a <- ggplot(b, aes(x=GuidePosition, y=b[[g]])) +
    geom_point(color = "darkgray") +
    geom_point(data=h, inherit.aes = F, 
      aes(x=GuidePosition,y=h[[g]]),color="red") +
    geom_hline(yintercept= 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept= -1, linetype = "dashed", color = "black") +
    #geom_text_repel(data=subset(
    #  b, g < (c-f) & guide_ctrl == "Screen"), #| log2fc_adj > f),
    #  aes(GuidePosition,g,label=guide),
    #  vjust=0,hjust=0,size=3, max.overlaps = 20) +
    #add enhancer rectangles with low alpha
    geom_rect(data = enhancer_key, inherit.aes = F, aes(
      xmin = Start, xmax = End, ymin = -Inf,ymax = Inf), alpha = 0.3, fill = "lightgray") +
    #theme_bw() + labs(x="Chr.9", y="Guide z-score") +
    theme_bw() + labs(x="Chr.9", y="Guide z-score") +
    theme(text = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    xlim(22064500,22163000)
  return(a)
}


pltplot <- function(a,b,g) {
  y <- b %>% 
    summarize(mean=mean(g),
              min=min(g),
              max= max(g),
              sd = sd(g))
  c=as.numeric(y[1])
  d = as.numeric(y[2])
  e = as.numeric(y[3])
  f = as.numeric(y[4])
  a <- ggplot(b, aes(x=GuidePosition, y=g, color=rep_id)) +
    geom_point() +
    # geom_point(data=subset(
    #   b, log2fc_adj < d),
    #   color="red") +
    # geom_point(data=subset(
    #   b, log2fc_adj > e),
    #   color="blue") +
    #geom_hline(yintercept=d, linetype = "dashed", color = "black") +
    #geom_hline(yintercept=e, linetype = "dashed", color = "black") +
    geom_hline(yintercept=1, linetype = "dashed", color = "black") +
    geom_hline(yintercept=-1, linetype = "dashed", color = "black") +
    #geom_text_repel(data=subset(
    #  b, log2fc_adj < d | log2fc_adj > e),
    #  aes(GuidePosition,log2fc_adj,label=guide),
    #  vjust=0,hjust=0,size=2) +
    #add enhancer rectangles with low alpha
    geom_rect(data = enhancer_key, inherit.aes = F, aes(
      xmin = Start, xmax = End, ymin = -Inf,ymax = Inf), alpha = 0.3, fill = "lightgray") +
    theme_bw() + labs(x="Chr.9", y="Guide z-score", color = "Plate") +
    theme(text = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    xlim(22064500,22163000)
}

####QC plots####

test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2B.mean) + facet_wrap(~plate) + ggtitle("CDKN2B")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_cdkn2b.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2B.mean) + ggtitle("CDKN2B")
ggsave(filename= file.path(save_path_batch, "smc_pooled_cdkn2b.png"), plot=test3, device=png, width = 10, height = 8)

test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2A.mean) + facet_wrap(~plate) + ggtitle("CDKN2A")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_cdkn2a.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2A.mean) + ggtitle("CDKN2A")
ggsave(filename= file.path(save_path_batch, "smc_pooled_cdkn2a.png"), plot=test3, device=png, width = 10, height = 8)

test3 <- pltplot(test3, smc_merged,smc_merged$DMRTA1.mean) + facet_wrap(~plate) + ggtitle("DMRTA1")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_dmrta1.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$DMRTA1.mean) + ggtitle("DMRTA1")
ggsave(filename= file.path(save_path_batch, "smc_pooled_dmrta1.png"), plot=test3, device=png, width = 10, height = 8)

test3 <- pltplot(test3, smc_merged,smc_merged$MTAP.mean) + facet_wrap(~plate) + ggtitle("MTAP")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_mtap.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$MTAP.mean) + ggtitle("MTAP")
ggsave(filename= file.path(save_path_batch, "smc_pooled_mtap.png"), plot=test3, device=png, width = 10, height = 8)

test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2B.AS1.mean) + facet_wrap(~plate) + ggtitle("CDKN2B-AS1")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_anril.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$CDKN2B.AS1.mean) + ggtitle("CDKN2B-AS1")
ggsave(filename= file.path(save_path_batch, "smc_pooled_anril.png"), plot=test3, device=png, width = 10, height = 8)

test3 <- pltplot(test3, smc_merged,smc_merged$MIR31HG.mean) + facet_wrap(~plate) + ggtitle("MIR31HG")
test3
ggsave(filename= file.path(save_path_batch, "smc_batch_mir31hg.png"), plot=test3, device=png, width = 10, height = 8)
test3 <- pltplot(test3, smc_merged,smc_merged$MIR31HG.mean) + ggtitle("MIR31HG")
ggsave(filename= file.path(save_path_batch, "smc_pooled_mir31hg.png"), plot=test3, device=png, width = 10, height = 8)


####Locus track plots####
enhancer_plot <- enhancer_plot + loc_lims
snp_plot <- snp_plot + loc_lims

cdkn2a_plot <- nnfplot(cdkn2a_plot,smc_merged,"CDKN2A.mean")
z <- cdkn2a_plot + ggtitle("CDKN2A") + loc_lims 
q <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#z <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_cdkn2a_locus_labs.png"), plot=q, device=png, width = 10, height = 8)


cdkn2b_plot <- nnfplot(cdkn2b_plot,smc_merged,"CDKN2B.mean")
y <- cdkn2b_plot + ggtitle("CDKN2B") + loc_lims
q <- y/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#y <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_cdkn2b_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

anril_plot <- nnfplot(anril_plot,smc_merged,"CDKN2B.AS1.mean")
x <- anril_plot + ggtitle("CDKN2B-AS1") + loc_lims
q <- x/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#x <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_anril_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

dmrta1_plot <- nnfplot(dmrta1_plot,smc_merged,"DMRTA1.mean")
w <- dmrta1_plot + ggtitle("DMRTA1") + loc_lims
q <- w/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#w <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_dmrta1_locus_labs.png"), plot=q, device=png, width = 10, height = 8)


mtap_plot <- nnfplot(mtap_plot,smc_merged,"MTAP.mean")
v <- mtap_plot +ggtitle("MTAP") + loc_lims
q <- v/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#v <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_mtap_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

mir31hg_plot <- nnfplot(mir31hg_plot,smc_merged,"MIR31HG.mean")
t <- mir31hg_plot +ggtitle("MIR31HG") + loc_lims
q <- t/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#v <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_figs, "smc_mir31hg_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

#plot all tracks together
u <- z/y/x/w/v/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=7, heights = c(1,1,1,1, 1, 0.3,0.3))
ggsave(filename= file.path(save_path_figs, "smc_1_all_tracks2.png"), plot=u, device=png, width = 24, height = 20)
ggsave(filename= file.path(save_path_figs, "smc_1_all_tracks2.svg"), plot=u, width = 24, height = 18)


####whole screen plots####

enhancer_plot <- enhancer_plot + screen_lims
snp_plot <- snp_plot + screen_lims

cdkn2a_plot <- nnfplot(cdkn2a_plot,smc_merged,"CDKN2A.mean")
z <- cdkn2a_plot + ggtitle("CDKN2A") + screen_lims 
q <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#z <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_cdkn2a_locus_labs.png"), plot=q, device=png, width = 10, height = 8)


cdkn2b_plot <- nnfplot(cdkn2b_plot,smc_merged,"CDKN2B.mean")
y <- cdkn2b_plot + ggtitle("CDKN2B") + screen_lims
q <- y/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#y <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_cdkn2b_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

anril_plot <- nnfplot(anril_plot,smc_merged,"CDKN2B.AS1.mean")
x <- anril_plot + ggtitle("CDKN2B-AS1") + screen_lims
q <- x/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#x <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_anril_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

dmrta1_plot <- nnfplot(dmrta1_plot,smc_merged,"DMRTA1.mean")
w <- dmrta1_plot + ggtitle("DMRTA1") + screen_lims
q <- w/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#w <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_dmrta1_locus_labs.png"), plot=q, device=png, width = 10, height = 8)


mtap_plot <- nnfplot(mtap_plot,smc_merged,"MTAP.mean")
v <- mtap_plot +ggtitle("MTAP") + screen_lims
q <- v/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#v <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_mtap_locus_labs.png"), plot=q, device=png, width = 10, height = 8)

mir31hg_plot <- nnfplot(mir31hg_plot,smc_merged,"MIR31HG.mean")
t <- mir31hg_plot +ggtitle("MIR31HG") + screen_lims
q <- t/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
q
#v <- z/enhancer_plot/snp_plot + plot_layout(ncol = 1, heights = c(1,0.2,0.2))
ggsave(filename= file.path(save_path_screen, "smc_mir31hg_locus_labs.png"), plot=q, device=png, width = 10, height = 8)


#plot all tracks together
u <- z/y/x/w/v/enhancer_plot/snp_plot + plot_layout(ncol = 1, nrow=7, heights = c(1,1,1,1, 1, 0.3,0.3))
ggsave(filename= file.path(save_path_screen, "smc_1_all_tracks2.png"), plot=u, device=png, width = 24, height = 20)
ggsave(filename= file.path(save_path_screen, "smc_1_all_tracks2.svg"), plot=u, width = 24, height = 18)








