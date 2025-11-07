#Script for z-score normalization of 9p21 potential target genes
#plots QC for read depth and gene counts

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
save_path_figs = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/locus_plots"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/datasets"
save_path_deseq = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/input"
save_path_batch = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/batch_effect"
save_path_batch_filtered = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/batch_effect/control_removed"

screen_lims = xlim(22064500,22163000)

loc_lims = xlim(22064500,22132000)

crop_lims = xlim(22095000,22120000)

#read in smc data by replicate

dxl04_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl04/dxl04_smc_rep1_tpm.csv")
dxl04_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl04/dxl04_smc_rep2_tpm.csv")

dxl05_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl05/dxl05_smc_rep1_tpm.csv")
dxl05_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl05/dxl05_smc_rep2_tpm.csv")

dxl02_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl02/dxl02_smc_rep1_tpm.csv")
dxl02_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl02/dxl02_smc_rep2_tpm.csv")

dxl03_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl03/dxl03_smc_rep1_tpm.csv")
dxl03_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl03/dxl03_smc_rep2_tpm.csv")

dxl01_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl01/dxl01_smc_rep1_tpm.csv")
dxl01_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxl01/dxl01_smc_rep2_tpm.csv")

dxk99_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk99/dxk99_smc_rep1_tpm.csv")
dxk99_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk99/dxk99_smc_rep2_tpm.csv")

dxk98_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk98/dxk98_smc_rep1_tpm.csv")
dxk98_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk98/dxk98_smc_rep2_tpm.csv")

dxk97_smc1 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk97/dxk97_smc_rep1_tpm.csv")
dxk97_smc2 <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/plate_qcs_and_reads/smc_dxk97/dxk97_smc_rep2_tpm.csv")

###normalize counts of each gene to the mean and standard deviation of the NT control###

means <- dxl04_smc1 %>% filter(is.na(depth)==F & is.na(guide_id)==F) %>% summarise(
  anril_mean = mean(CDKN2B.AS1),
  anril_sd = sd(CDKN2B.AS1),
  mtap_mean = mean(MTAP),
  mtap_sd = sd(MTAP),
  cdkn2a_mean = mean(CDKN2A),
  cdkn2a_sd = sd(CDKN2A),
  cdkn2b_mean = mean(CDKN2B),
  cdkn2b_sd = sd(CDKN2B),
  dmrta1_mean = mean(DMRTA1),
  dmrta1_sd = sd(DMRTA1),
  mir31hg_mean = mean(MIR31HG),
  mir31hg_sd = sd(MIR31HG)
  
)

#a=input data
#b=output file with normalized values
#perform z-score normalization of gene measurements for each sample
cpm_normalize <- function(a) {
  means <- a %>% filter(is.na(depth)==F & is.na(guide_id)==F) %>% summarise(
    anril_mean = mean(CDKN2B.AS1),
    anril_sd = sd(CDKN2B.AS1),
    mtap_mean = mean(MTAP),
    mtap_sd = sd(MTAP),
    cdkn2a_mean = mean(CDKN2A),
    cdkn2a_sd = sd(CDKN2A),
    cdkn2b_mean = mean(CDKN2B),
    cdkn2b_sd = sd(CDKN2B),
    dmrta1_mean = mean(DMRTA1),
    dmrta1_sd = sd(DMRTA1),
    mir31hg_mean = mean(MIR31HG),
    mir31hg_sd = sd(MIR31HG)
  )
  b <- a
  b <- b %>% mutate(
    CDKN2B.AS1 = (CDKN2B.AS1-means$anril_mean)/means$anril_sd,
    MTAP = (MTAP-means$mtap_mean)/means$mtap_sd,
    CDKN2A = (CDKN2A-means$cdkn2a_mean)/means$cdkn2a_sd,
    CDKN2B = (CDKN2B-means$cdkn2b_mean)/means$cdkn2b_sd,
    DMRTA1 = (DMRTA1-means$dmrta1_mean)/means$dmrta1_sd,
    MIR31HG = (MIR31HG-means$mir31hg_mean)/means$mir31hg_sd
  )
  return(b)
}

dxl04_smc1 <- cpm_normalize(dxl04_smc1)
dxl04_smc2 <- cpm_normalize(dxl04_smc2)
dxl05_smc1 <- cpm_normalize(dxl05_smc1)
dxl05_smc2 <- cpm_normalize(dxl05_smc2)
dxl03_smc1 <- cpm_normalize(dxl03_smc1)
dxl03_smc2 <- cpm_normalize(dxl03_smc2)
dxl02_smc1 <- cpm_normalize(dxl02_smc1)
dxl02_smc2 <- cpm_normalize(dxl02_smc2)
dxl01_smc1 <- cpm_normalize(dxl01_smc1)
dxl01_smc2 <- cpm_normalize(dxl01_smc2)
dxk99_smc1 <- cpm_normalize(dxk99_smc1)
dxk99_smc2 <- cpm_normalize(dxk99_smc2)
dxk98_smc1 <- cpm_normalize(dxk98_smc1)
dxk98_smc2 <- cpm_normalize(dxk98_smc2)
dxk97_smc1 <- cpm_normalize(dxk97_smc1)
dxk97_smc2 <- cpm_normalize(dxk97_smc2)

####assign guide names for use with DESEQ2 later on VERY IMPORTANT####
smc_reps <- rbind(dxl04_smc1,dxl05_smc1,dxl02_smc1,dxl03_smc1,
                  dxl01_smc1,dxk99_smc1,dxk98_smc1, dxk97_smc1)

smc_reps <- smc_reps %>% mutate(
  guide_name = case_when(guide == "NEG_CONTROL" ~ paste0(guide_id, "_rep", rep_id),
                         sequence != "" ~ paste0(guide_id, "_rep1"))
)

smc_reps2 <- rbind(dxl04_smc2,dxl05_smc2,dxl02_smc2,dxl03_smc2,
                   dxl01_smc2,dxk99_smc2, dxk98_smc2, dxk97_smc2)

smc_reps2 <- smc_reps2 %>% mutate(
  guide_name = case_when(guide == "NEG_CONTROL" ~ paste0(guide_id, "_rep", rep_id),
                         sequence != "" ~ paste0(guide_id, "_rep2"))
)

smc_merged <- rbind(smc_reps,smc_reps2)


#smc_merged <- rbind(dxl04_smc1,dxl04_smc2,dxl05_smc1,dxl05_smc2,dxl02_smc1,dxl02_smc2,dxl03_smc1,dxl03_smc2,
#                    dxl01_smc1,dxl01_smc2,dxk99_smc1,dxk99_smc2,dxk98_smc1,dxk98_smc2, dxk97_smc1, dxk97_smc2)

#remove random "sample number" column and move guide rep column to front
smc_merged <- smc_merged[,-1] %>% relocate(guide_name)

#save file of all readcounts and metadata for regression
write.csv(smc_merged, paste0(save_path_files,"/plate_normalized_cpm_reads_all_input_all_ctrl.csv"))

smc_merged <- smc_merged %>% filter(!is.na(guide_id))
smc_merged <- smc_merged %>% filter(!is.na(depth))
write.csv(smc_merged, paste0(save_path_deseq,"/plate_normalized_cpm_reads_all_input.csv"))

####create DF with both replicates for each guide in the same row, to compare guide agreement####

smc_reps <- merge(smc_reps,smc_reps2,by=c("guide_id","plate"))
smc_reps <- smc_reps %>% filter(!is.na(guide_id))

rep_plot <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps, aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(aes(color = guide_ctrl.x), alpha=0.6) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    #(aes(color=guide_ctrl.x)
    geom_smooth(method = lm, color="grey") + theme_bw() +
    scale_color_manual(values = c("blue","skyblue","gray20")) +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}

nongene_rep_plot <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps, aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(aes(color=plate), alpha=0.6) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    #geom_smooth(method = lm, color="grey") + 
    theme_bw() +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}


nongene_rep_plot_2 <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps, aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(aes(color=guide_ctrl.x), alpha=0.6) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    geom_smooth(method = lm, color="grey") + theme_bw() +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}


#lm_eq <- function(z) {substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
lm_eq <- function(z) {substitute(italic(r)^2~"="~r2,
                                 
                                 list(        #a = as.character(as.data.frame(format(coef(z)[1], digits = 3))),
                                              
                                              #b = as.character(as.data.frame(format(coef(z)[2], digits = 3))),
                                              
                                              r2 = format(summary(z)$r.squared, digits = 3)))
}


####Plot correlation of 9p21 genes across replicates as well as all other features####

model_anril <-lm(`CDKN2B.AS1.y`~`CDKN2B.AS1.x`,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_mtap <-lm(MTAP.y~MTAP.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_cdkn2b <-lm(CDKN2B.y~CDKN2B.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_cdkn2a <-lm(CDKN2A.y~CDKN2A.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_dmrta1 <-lm(DMRTA1.y~DMRTA1.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_mir31hg <-lm(MIR31HG.y~MIR31HG.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

eqs <- list(
  eq_anril = lm_eq(model_anril),
  eq_mtap = lm_eq(model_mtap),
  eq_cdkn2b = lm_eq(model_cdkn2b),
  eq_cdkn2a = lm_eq(model_cdkn2a),
  eq_dmrta1 = lm_eq(model_dmrta1),
  eq_mir31hg = lm_eq(model_mir31hg)
)

plot_save <- rep_plot("CDKN2B.AS1.x", "CDKN2B.AS1.y") +
  labs(x="CDKN2B-AS1 z-score in first replicate", y="CDKN2B-AS1 z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_anril))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_anril.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.AS1.x", "CDKN2B.AS1.y") + facet_wrap(~plate) +
  labs(x="CDKN2B-AS1 z-score in first replicate", y="CDKN2B-AS1 z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_anril_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.x", "CDKN2B.y") +
  labs(x="CDKN2B z-score in first replicate", y="CDKN2B z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_cdkn2b))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_cdkn2b.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.x", "CDKN2B.y") + facet_wrap(~plate) +
  labs(x="CDKN2B z-score in first replicate", y="CDKN2B z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_cdkn2b_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2A.x", "CDKN2A.y") +
  labs(x="CDKN2A z-score in first replicate", y="CDKN2A z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_cdkn2a))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_cdkn2a.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2A.x", "CDKN2A.y") + facet_wrap(~plate) +
  labs(x="CDKN2A z-score in first replicate", y="CDKN2A z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_cdkn2a_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MTAP.x", "MTAP.y") +
  labs(x="MTAP z-score in first replicate", y="MTAP z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_mtap))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_mtap.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MTAP.x", "MTAP.y") + facet_wrap(~plate) +
  labs(x="MTAP z-score in first replicate", y="MTAP z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_mtap_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("DMRTA1.x", "DMRTA1.y") +
  labs(x="DMRTA1 z-score in first replicate", y="DMRTA1 z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_dmrta1))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_dmrta1.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("DMRTA1.x", "DMRTA1.y") + facet_wrap(~plate) +
  labs(x="DMRTA1 z-score in first replicate", y="DMRTA1 z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_dmrta1_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MIR31HG.x", "MIR31HG.y") +
  labs(x="MIR31HG z-score in first replicate", y="MIR31HG z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_mir31hg))
plot_save
ggsave(filename= file.path(save_path_batch, "smc_replicates_mir31hg.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MIR31HG.x", "MIR31HG.y") + facet_wrap(~plate) +
  labs(x="MIR31HG z-score in first replicate", y="MIR31HG z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_mir31hg_plates.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot("depth.x", "depth.y") +
  labs(x="Total reads in first replicate", y="Total reads in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_read_depth.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot_2("depth.x", "depth.y") + facet_wrap(~plate) +
  labs(x="Total reads in first replicate", y="Total reads in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_read_depth_plate.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot("mtper.x", "mtper.y") +
  labs(x="Mito fraction in first replicate", y="Mito fraction in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_mito_frac.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot("riboper.x", "riboper.y") +
  labs(x="Ribosomal fraction in first replicate", y="Ribosomal fraction in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch, "smc_replicates_ribo_frac.png"), plot=plot_save, device=png, width = 10, height = 8)

#plot_save <- rep_plot("CDKN2B.AS1.x", "CDKN2B.AS1.y") +
#  labs(x="CDKN2B-AS1 z-score in first replicate", y="CDKN2B-AS1 z-score in second replicate", color = "Guide target")
#ggsave(filename= file.path(save_path_batch, "smc_replicates_anril.png"), plot=plot_save, device=png, width = 10, height = 8)


####Plot data again but remove controls####

smc_reps <- smc_reps %>% filter(guide_ctrl.x == "Screen")

model_anril <-lm(`CDKN2B.AS1.y`~`CDKN2B.AS1.x`,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_mtap <-lm(MTAP.y~MTAP.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_cdkn2b <-lm(CDKN2B.y~CDKN2B.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_cdkn2a <-lm(CDKN2A.y~CDKN2A.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_dmrta1 <-lm(DMRTA1.y~DMRTA1.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

model_mir31hg <-lm(MIR31HG.y~MIR31HG.x,data=smc_reps %>% filter(guide_ctrl.x == "Screen"))

eqs <- list(
  eq_anril = lm_eq(model_anril),
  eq_mtap = lm_eq(model_mtap),
  eq_cdkn2b = lm_eq(model_cdkn2b),
  eq_cdkn2a = lm_eq(model_cdkn2a),
  eq_dmrta1 = lm_eq(model_dmrta1),
  eq_mir31hg = lm_eq(model_mir31hg)
)

rep_plot <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps %>% filter(guide_ctrl.x == "Screen"), aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(color = "gray20", alpha=0.6) + coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    #(aes(color=guide_ctrl.x)
    geom_smooth(method = lm, color="grey") + theme_bw() +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}

nongene_rep_plot <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps %>% filter(guide_ctrl.x == "Screen"), aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(aes(color=plate), alpha=0.6) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    #geom_smooth(method = lm, color="grey") + 
    theme_bw() +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}


nongene_rep_plot_2 <- function(xaxis, yaxis) {
  zscore_plot <- ggplot(smc_reps %>% filter(guide_ctrl.x == "Screen"), aes(x=.data[[xaxis]], y=.data[[yaxis]])) +
    geom_point(aes(color=guide_ctrl.x), alpha=0.6) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    geom_smooth(method = lm, color="grey") + theme_bw() +
    theme(text = element_text(size = 15),
          title = element_text(size = 20),
          legend.text = element_text(size=15),
          axis.title.x=element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  return(zscore_plot)
}

plot_save <- rep_plot("CDKN2B.AS1.x", "CDKN2B.AS1.y") +
  labs(x="CDKN2B-AS1 z-score in first replicate", y="CDKN2B-AS1 z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_anril))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_anril_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.AS1.x", "CDKN2B.AS1.y") + facet_wrap(~plate) +
  labs(x="CDKN2B-AS1 z-score in first replicate", y="CDKN2B-AS1 z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_anril_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.x", "CDKN2B.y") +
  labs(x="CDKN2B z-score in first replicate", y="CDKN2B z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_cdkn2b))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_cdkn2b_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2B.x", "CDKN2B.y") + facet_wrap(~plate) +
  labs(x="CDKN2B z-score in first replicate", y="CDKN2B z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_cdkn2b_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2A.x", "CDKN2A.y") +
  labs(x="CDKN2A z-score in first replicate", y="CDKN2A z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_cdkn2a))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_cdkn2a_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("CDKN2A.x", "CDKN2A.y") + facet_wrap(~plate) +
  labs(x="CDKN2A z-score in first replicate", y="CDKN2A z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_cdkn2a_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MTAP.x", "MTAP.y") +
  labs(x="MTAP z-score in first replicate", y="MTAP z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_mtap))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_mtap_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MTAP.x", "MTAP.y") + facet_wrap(~plate) +
  labs(x="MTAP z-score in first replicate", y="MTAP z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_mtap_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("DMRTA1.x", "DMRTA1.y") +
  labs(x="DMRTA1 z-score in first replicate", y="DMRTA1 z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_dmrta1))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_dmrta1_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("DMRTA1.x", "DMRTA1.y") + facet_wrap(~plate) +
  labs(x="DMRTA1 z-score in first replicate", y="DMRTA1 z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_dmrta1_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MIR31HG.x", "MIR31HG.y") +
  labs(x="MIR31HG z-score in first replicate", y="MIR31HG z-score in second replicate", color = "Guide target") +
  annotate(geom="text", x=1.5, y=-3, label = as.call(eqs$eq_mir31hg))
plot_save
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_mir31hg_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- rep_plot("MIR31HG.x", "MIR31HG.y") + facet_wrap(~plate) +
  labs(x="MIR31HG z-score in first replicate", y="MIR31HG z-score in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_mir31hg_plates_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)


plot_save <- nongene_rep_plot("depth.x", "depth.y") +
  labs(x="Total reads in first replicate", y="Total reads in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_read_depth_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot_2("depth.x", "depth.y") + facet_wrap(~plate) +
  labs(x="Total reads in first replicate", y="Total reads in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_read_depth_plate_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot("mtper.x", "mtper.y") +
  labs(x="Mito fraction in first replicate", y="Mito fraction in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_mito_frac_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)

plot_save <- nongene_rep_plot("riboper.x", "riboper.y") +
  labs(x="Ribosomal fraction in first replicate", y="Ribosomal fraction in second replicate", color = "Guide target")
ggsave(filename= file.path(save_path_batch_filtered, "smc_replicates_ribo_frac_filtered.png"), plot=plot_save, device=png, width = 10, height = 8)


smc_reps_avg <- smc_reps %>% group_by(plate) %>% summarize(
  mito_frac_x = mean(mtper.x, na.rm=T),
  mito_frac_y = mean(mtper.y, na.rm=T),
  mito_sd_x = sd(mtper.x, na.rm=T),
  mito_sd_y = sd(mtper.y, na.rm=T),
  minx = mito_frac_x-mito_sd_x,
  maxx = mito_frac_x+mito_sd_x,
  miny = mito_frac_y-mito_sd_y,
  maxy = mito_frac_y+mito_sd_y,
  ribo_frac_x = mean(riboper.x, na.rm=T),
  ribo_frac_y = mean(riboper.y, na.rm=T),
  ribo_sd_x = sd(riboper.x, na.rm=T),
  ribo_sd_y = sd(riboper.y, na.rm=T),
  rminx = ribo_frac_x-ribo_sd_x,
  rmaxx = ribo_frac_x+ribo_sd_x,
  rminy = ribo_frac_y-ribo_sd_y,
  rmaxy = ribo_frac_y+ribo_sd_y
)

zscore_plot <- ggplot(smc_reps_avg, aes(x=mito_frac_x, y=mito_frac_y)) +
  geom_errorbar(aes(ymin = miny,ymax = maxy)) + 
  geom_errorbarh(aes(xmin = minx,xmax = maxx)) + 
  geom_point(aes(color=plate), size=3) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size=15),
        axis.title.x=element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(x="Mitochondrial fraction in replicate 1", y="Mitochondrial fraction in replicate 2", color = "Virus plate")
zscore_plot
ggsave(filename= file.path(save_path_batch, "smc_replicates_mito_sd.png"), plot=zscore_plot, device=png, width = 8, height = 6)


zscore_plot <- ggplot(smc_reps_avg, aes(x=ribo_frac_x, y=ribo_frac_y)) +
  geom_errorbar(aes(ymin = rminy,ymax = rmaxy)) + 
  geom_errorbarh(aes(xmin = rminx,xmax = rmaxx)) + 
  geom_point(aes(color=plate), size=3) + #coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_bw() +
  theme(text = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size=15),
        axis.title.x=element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15)) +
  labs(x="Ribosomal fraction in replicate 1", y="Ribosomal fraction in replicate 2", color = "Virus plate")
zscore_plot
ggsave(filename= file.path(save_path_batch, "smc_replicates_ribo_sd.png"), plot=zscore_plot, device=png, width = 8, height = 6)



