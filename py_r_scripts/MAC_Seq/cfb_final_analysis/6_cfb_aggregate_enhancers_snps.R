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
save_path_figs = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/cfb/zscore_normalized/locus_plots"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/cfb/zscore_normalized/datasets"
save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/cfb/zscore_normalized/sliding_window"
save_path_snps = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/cfb/zscore_normalized/enhancer_and_snp_plots"
save_path_tss = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/regression_output/cfb/zscore_normalized/tss"

export_fig2 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/screen_fig"
export_fig4 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/supp_fig4"
export_fig2 = "/Volumes/broad_mcl/members_dir/bschmand/9p21_manuscript_figure_data/fig2"

cfb_merged <- read.csv(paste0(save_path_files,"/plate_normalized_guide_means_output.csv"))
#merged_sample_data <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/deseq2/datasets/plate_normalized_guide_means_output.csv")
cfb_merged <- cfb_merged[,-1]

####make enhancer and SNP plots for use with other graphics####
#enhancer_key <- read.csv("raw_data/enhancer_key.csv")
enhancer_key <- read.csv(paste0(input_path,"/enhancer_key_controls.csv"))
#enhancer_key_only <- read.csv("raw_data/enhancer_key.csv")
enhancer_key_only <- read.csv(paste0(input_path,"/enhancer_key.csv"))
enhancer_key$Start <- as.numeric(gsub(",","",enhancer_key$Start))
enhancer_key$End <-  as.numeric(gsub(",","",enhancer_key$End))
enhancer_key$Mid <-  as.numeric(gsub(",","",enhancer_key$Mid))

#snp_key <- read.csv("raw_data/SNP_key.csv",fileEncoding="latin1")
snp_key <- read.csv(paste0(input_path,"/SNP_key.csv"),fileEncoding="latin1")
snp_key$loc <- as.numeric(snp_key$loc)

####sort data and plot violin by enhancer####

#plot only guides labeled by enhancer

test_data <- cfb_merged %>% filter(GuidePosition < 22132000) %>%
  filter(grepl("E|Neg|NEG", guide))
test <- ggplot(test_data, aes(x=guide, y=DMRTA1.mean))
test <-  test + geom_violin() + geom_point(aes(color=plate)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
test

test_data <- cfb_merged %>% filter(GuidePosition < 22132000) %>%
  filter(!grepl("rs", guide))

test <- ggplot(test_data, aes(x=guide, y=DMRTA1))
test <-  test + geom_violin() + geom_point()
test

#select all guides within an enhancer range

dmrta1_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "DMRTA1") ##%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  dmrta1_enhancers <- rbind(dmrta1_enhancers,y)
}

means <- dmrta1_enhancers %>% group_by(near) %>% 
  summarise(fc_mean=mean(DMRTA1.mean),
            fc_min = min(DMRTA1.mean),
            fc_max = max(DMRTA1.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(fc_mean = ifelse(near == "loc.", 20,fc_mean)
  )
dmrta1_enhancers <- merge(dmrta1_enhancers, means)
z <- subset(dmrta1_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(DMRTA1),mean=mean(DMRTA1))
x <- as.numeric(z[1])
#y <- as.numeric(z[2])

test_data <- dmrta1_enhancers %>% mutate(near = fct_reorder(near, fc_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476 | GuidePosition > 22447500)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=DMRTA1.mean))
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=fc_min, ymax=fc_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "DMRTA1 gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "dmrta1_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)

test <- dmrta1_enhancers[5:722,1:32]
test <- test[!duplicated(test$guide_id),]
write.csv(dmrta1_enhancers[5:722,1:32], paste0(save_path_snps,"/enhancer_effects_for_prism.csv"))


#repeat for mtap
mtap_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "MTAP") ##%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  mtap_enhancers <- rbind(mtap_enhancers,y)
}

means <- mtap_enhancers %>% group_by(near) %>% 
  summarise(fc_mean=mean(MTAP.mean),
            fc_min = min(MTAP.mean),
            fc_max = max(MTAP.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(fc_mean = ifelse(near == "loc.", 20,fc_mean)
  )
mtap_enhancers <- merge(mtap_enhancers, means)
z <- subset(mtap_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(MTAP.mean),mean=mean(MTAP.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- mtap_enhancers %>% mutate(near = fct_reorder(near, fc_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=MTAP.mean)) + ggtitle("MTAP gene expression by enhancer")
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=fc_min, ymax=fc_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "MTAP gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "mtap_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)

#export unique enhancer-snp intersections
enhancers_in_snps <- mtap_enhancers %>% select(near,guide) %>% filter(grepl("rs",guide))
enhancers_in_snps <- enhancers_in_snps[!duplicated(enhancers_in_snps$guide),]
write.csv(enhancers_in_snps, paste0(save_path_files,"/enhancer_snp_intersection.csv"))


#for cdkn2a
cdkn2a_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "CDKN2A") #%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  cdkn2a_enhancers <- rbind(cdkn2a_enhancers,y)
}

means <- cdkn2a_enhancers %>% group_by(near) %>% 
  summarise(snp_mean=mean(CDKN2A.mean),
            snp_min = min(CDKN2A.mean),
            snp_max = max(CDKN2A.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean)
  )
cdkn2a_enhancers <- merge(cdkn2a_enhancers, means)
z <- subset(cdkn2a_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2A.mean),mean=mean(CDKN2A.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- cdkn2a_enhancers %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=CDKN2A.mean)) + ggtitle("CDKN2A gene expression by enhancer")
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "CDKN2A gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "cdkn2a_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)

#for cdkn2b
cdkn2b_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "CDKN2B") #%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  cdkn2b_enhancers <- rbind(cdkn2b_enhancers,y)
}

means <- cdkn2b_enhancers %>% group_by(near) %>% 
  summarise(fc_mean=mean(CDKN2B.mean),
            fc_min = min(CDKN2B.mean),
            fc_max = max(CDKN2B.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(fc_mean = ifelse(near == "loc.", 20,fc_mean)
  )
cdkn2b_enhancers <- merge(cdkn2b_enhancers, means)
z <- subset(cdkn2b_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2B.mean),mean=mean(CDKN2B.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- cdkn2b_enhancers %>% mutate(near = fct_reorder(near, fc_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=CDKN2B.mean)) + ggtitle("CDKN2B gene expression by enhancer")
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=fc_min, ymax=fc_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "CDKN2B gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "cdkn2b_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)

#for anril
anril_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "CDKN2B-AS1") #%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  anril_enhancers <- rbind(anril_enhancers,y)
}

means <- anril_enhancers %>% group_by(near) %>% 
  summarise(fc_mean=mean(CDKN2B.AS1.mean),
            fc_min = min(CDKN2B.AS1.mean),
            fc_max = max(CDKN2B.AS1.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(fc_mean = ifelse(near == "loc.", 20,fc_mean)
  )
anril_enhancers <- merge(anril_enhancers, means)
z <- subset(anril_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2B.AS1.mean),mean=mean(CDKN2B.AS1.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- anril_enhancers %>% mutate(near = fct_reorder(near, fc_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=CDKN2B.AS1.mean))
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=fc_min, ymax=fc_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "ANRIL gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "anril_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)

#for mir31hg
mir31hg_enhancers <- cfb_merged %>% filter(guide == "Neg-9p21 regions" | guide == "MIR31HG") #%>% mutate(near = guide)

for (i in 1:nrow(enhancer_key_only)) {
  z <- enhancer_key_only[i,2]
  x <- enhancer_key_only[i,3]
  y <- cfb_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  mir31hg_enhancers <- rbind(mir31hg_enhancers,y)
}

means <- mir31hg_enhancers %>% group_by(near) %>% 
  summarise(fc_mean=mean(MIR31HG.mean),
            fc_min = min(MIR31HG.mean),
            fc_max = max(MIR31HG.mean)) %>%
  #assign neg 9p21 regions an arbitrary "high mean" to plot it first
  mutate(fc_mean = ifelse(near == "loc.", 20,fc_mean)
  )
mir31hg_enhancers <- merge(mir31hg_enhancers, means)
z <- subset(mir31hg_enhancers, guide=="Neg-9p21 regions") %>% summarize(sd = sd(MIR31HG.mean),mean=mean(MIR31HG.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- mir31hg_enhancers %>% mutate(near = fct_reorder(near, fc_mean, .desc=T)) %>% 
  filter(GuidePosition < 22136476)
#test_data <- test_data[-1,]
enhancer_plot <- ggplot(test_data, aes(x=near, y=MIR31HG.mean))
enhancer_plot <-  enhancer_plot + geom_violin() + #geom_pointrange(aes(ymin=fc_min, ymax=fc_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "MIR31HG gene expression by enhancer", y="Guide z-score", x="Targeted region")
enhancer_plot

ggsave(filename= file.path(save_path_snps, "mir31hg_split_enhancers.png"), 
       plot=enhancer_plot, device=png, width = 10, height = 5)


####plot SNPs and controls, aggregate 300bp####

#plot SNPs including all guides within 300bp window

#for dmrta1: group by SNPs
#start dataframe using the negative controls we will plot, 
#"near" will be the column name for identifying any nearby variant/feature

dmrta1_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  dmrta1_snps <- rbind(dmrta1_snps,y)
}
#group by SNP
means <- dmrta1_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(DMRTA1.mean),
            snp_min = min(DMRTA1.mean),
            snp_max = max(DMRTA1.mean))
dmrta1_snps <- merge(dmrta1_snps, means)

dmrta1_snps <- dmrta1_snps %>% arrange(near,snp_mean)
#write.csv(dmrta1_snps, paste0(save_path_files,"/dmrta1_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
dmrta1_snps <- dmrta1_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(dmrta1_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(DMRTA1.mean),mean=mean(DMRTA1.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- dmrta1_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000) #22132000
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=(DMRTA1.mean)))
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-x,x), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "DMRTA1 gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "dmrta1_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

#for cdkn2a: group by SNPs
cdkn2a_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  cdkn2a_snps <- rbind(cdkn2a_snps,y)
}
#group by SNP
means <- cdkn2a_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(CDKN2A.mean),
            snp_min = min(CDKN2A.mean),
            snp_max = max(CDKN2A.mean))
cdkn2a_snps <- merge(cdkn2a_snps, means)

cdkn2a_snps <- cdkn2a_snps %>% arrange(near,snp_mean)
#write.csv(cdkn2a_snps, paste0(save_path_files,"/cdkn2a_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
cdkn2a_snps <- cdkn2a_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(cdkn2a_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2A.mean),mean=mean(CDKN2A.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- cdkn2a_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000)
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=CDKN2A.mean))
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "CDKN2A gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "cdkn2a_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

#for anril: group by SNPs
anril_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  anril_snps <- rbind(anril_snps,y)
}
#group by SNP
means <- anril_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(CDKN2B.AS1.mean),
            snp_min = min(CDKN2B.AS1.mean),
            snp_max = max(CDKN2B.AS1.mean))
anril_snps <- merge(anril_snps, means)

anril_snps <- anril_snps %>% arrange(near,snp_mean)
#write.csv(anril_snps, paste0(save_path_files,"/anril_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
anril_snps <- anril_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(anril_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2B.AS1.mean),mean=mean(CDKN2B.AS1.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- anril_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000)
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=CDKN2B.AS1.mean))
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "ANRIL gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "anril_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

#for mtap: group by SNPs
mtap_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  mtap_snps <- rbind(mtap_snps,y)
}
#group by SNP
means <- mtap_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(MTAP.mean),
            snp_min = min(MTAP.mean),
            snp_max = max(MTAP.mean))
mtap_snps <- merge(mtap_snps, means)

mtap_snps <- mtap_snps %>% arrange(near,snp_mean)
#write.csv(mtap_snps, paste0(save_path_files,"/mtap_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
mtap_snps <- mtap_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(mtap_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(MTAP.mean),mean=mean(MTAP.mean))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- mtap_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000)
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=MTAP.mean))
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "MTAP gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "mtap_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

#for cdkn2b: group by SNPs
cdkn2b_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  cdkn2b_snps <- rbind(cdkn2b_snps,y)
}
#group by SNP
means <- cdkn2b_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(CDKN2B.mean),
            snp_min = min(CDKN2B.mean),
            snp_max = max(CDKN2B.mean))
cdkn2b_snps <- merge(cdkn2b_snps, means)

cdkn2b_snps <- cdkn2b_snps %>% arrange(near,snp_mean)
write.csv(cdkn2b_snps, paste0(save_path_files,"/cdkn2b_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
cdkn2b_snps <- cdkn2b_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(cdkn2b_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(CDKN2B),mean=mean(CDKN2B))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- cdkn2b_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000)
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=CDKN2B.mean)) #ggtitle("CDKN2B gene expression by variant")
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "CDKN2B gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "cdkn2b_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

#for mir31hg: group by SNPs
mir31hg_snps <- cfb_merged %>% filter(guide == "Neg-9p21 regions") #%>% mutate(near = guide)
for (i in 1:nrow(snp_key)) {
  z <- snp_key[i,2]
  y <- cfb_merged %>% filter((z-150)<GuidePosition & (z+150)>GuidePosition)
  y <- y %>% mutate(near = snp_key[i,1])
  mir31hg_snps <- rbind(mir31hg_snps,y)
}
#group by SNP
means <- mir31hg_snps %>% group_by(near) %>% 
  summarise(snp_mean=mean(MIR31HG.mean),
            snp_min = min(MIR31HG.mean),
            snp_max = max(MIR31HG.mean))
mir31hg_snps <- merge(mir31hg_snps, means)

mir31hg_snps <- mir31hg_snps %>% arrange(near,snp_mean)
write.csv(mir31hg_snps, paste0(save_path_files,"/mir31hg_snp_effects.csv"))

#assign neg 9p21 regions an arbitrary "high mean" to plot it first
mir31hg_snps <- mir31hg_snps %>% mutate(snp_mean = ifelse(near == "loc.", 20,snp_mean))
z <- subset(mir31hg_snps, guide=="Neg-9p21 regions") %>% summarize(sd = sd(MIR31HG),mean=mean(MIR31HG))
x <- as.numeric(z[1])
y <- as.numeric(z[2])

test_data <- mir31hg_snps %>% mutate(near = fct_reorder(near, snp_mean, .desc=T)) %>% 
  filter(GuidePosition < 22132000)
#test_data <- test_data[-1,]
snp_plot <- ggplot(test_data, aes(x=near, y=MIR31HG.mean)) #ggtitle("MIR31HG gene expression by variant")
snp_plot <-  snp_plot + geom_pointrange(aes(ymin=snp_min, ymax=snp_max)) +   
  geom_hline(yintercept=c(-1,1), linetype = "dashed", color = "black") +
  geom_point(aes(color=plate)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = "MIR31HG gene expression by variant", y="Guide z-score", x="Targeted region")
snp_plot

ggsave(filename= file.path(save_path_snps, "mir31hg_split_snps.png"), 
       plot=snp_plot, device=png, width = 10, height = 3)

####TSS plotting####


tss_select <- cfb_merged %>% filter(near %in% c("NT", "CDKN2B","CDKN2A","CDKN2B.AS1","DMRTA1","MTAP"))

#tss_means <- tss_select %>% group_by

tss_funct <- function(gene) {
  tss_input <- tss_select %>% filter(near %in% c("NT",gene))
  tss_input$near <- factor(tss_input$near, levels = c("NT", gene))
  t <- ggplot(tss_input, aes(x=near,y=.data[[paste0(gene, ".mean")]])) + geom_boxplot(aes(fill = near)) + geom_point(size = 2) +
    theme_classic() +
    theme(text = element_text(size = 20),
          #legend.text = element_blank(),
          legend.text = element_text(size=15),
          axis.title.x= element_blank(),
          axis.title.y=element_text(size = 15),
          axis.text.x = element_text(size = 15, angle = 70, hjust = 1),
          axis.text.y = element_text(size = 10)) +
    #remove legend
    guides(fill="none")
  return(t)
}

tss_plot <- tss_funct("MTAP") +
  labs(x="Treatment", y="Guide z-score") +
  scale_fill_manual(values = c("NT" = "gray90", "MTAP" = "lightskyblue2"), name = "Guide", labels = c("Scrambled control", "MTAP TSS")) +
  scale_x_discrete(labels = c("Scrambled\ncontrol", "MTAP TSS")) +
  stat_compare_means(label = "p.format", method = "t.test", label.x = 1.35, label.y = 1.3) +
  annotate("segment",x = 1.2, xend = 1.8, y = 1.2, yend = 1.2)
tss_plot

ggsave(filename= file.path(save_path_tss, "tss_mtap.png"), 
       plot=tss_plot, device=png, width = 2.5, height = 5)

ggsave(filename= file.path(export_fig4, "cfb_tss_mtap.svg"), 
       plot=tss_plot, width = 2.5, height = 5)

tss_plot <- tss_funct("DMRTA1") +
  labs(x="Treatment", y="Guide z-score")+
  scale_fill_manual(values = c("NT" = "gray90", "DMRTA1" = "lightskyblue2"), name = "Guide", labels = c("Scrambled control", "DMRTA1 TSS")) +
  scale_x_discrete(labels = c("Scrambled\ncontrol", "DMRTA1 TSS")) +
  stat_compare_means(label = "p.format", method = "t.test", label.x = 1.35, label.y = 0.45) +
  annotate("segment",x = 1.2, xend = 1.8, y = 0.4, yend = 0.4)
tss_plot

ggsave(filename= file.path(save_path_tss, "tss_dmrta1.png"), 
       plot=tss_plot, device=png, width = 2.5, height = 5)

ggsave(filename= file.path(export_fig4, "cfb_tss_dmrta1.svg"), 
       plot=tss_plot, width = 2.5, height = 5)

tss_plot <- tss_funct("CDKN2B") +
  labs(x="Treatment", y="Guide z-score") +
  scale_fill_manual(values = c("NT" = "gray90", "CDKN2B" = "lightskyblue2"), name = "Guide", labels = c("Scrambled control", "CDKN2B TSS")) +
  scale_x_discrete(labels = c("Scrambled\ncontrol", "CDKN2B TSS")) +
  stat_compare_means(label = "p.format", method = "t.test", label.x = 1.35, label.y = 1.3) +
  annotate("segment",x = 1.2, xend = 1.8, y = 1.2, yend = 1.2)
tss_plot

ggsave(filename= file.path(save_path_tss, "tss_cdkn2b.png"), 
       plot=tss_plot, device=png, width = 2.5, height = 5)

ggsave(filename= file.path(export_fig4, "cfb_tss_cdkn2b.svg"), 
       plot=tss_plot, width = 2.5, height = 5)


tss_plot <- tss_funct("CDKN2A") +
  labs(x="Treatment", y="Guide z-score")+
  scale_fill_manual(values = c("NT" = "gray90", "CDKN2A" = "lightskyblue2"), name = "Guide", labels = c("Scrambled control", "CDKN2A TSS")) +
  scale_x_discrete(labels = c("Scrambled\ncontrol", "CDKN2A TSS")) +
  stat_compare_means(label = "p.format", method = "t.test", label.x = 1.35, label.y = 1.1) +
  annotate("segment",x = 1.2, xend = 1.8, y = 1, yend = 1)
tss_plot

ggsave(filename= file.path(save_path_tss, "tss_cdkn2a.png"), 
       plot=tss_plot, device=png, width = 2.5, height = 5)

ggsave(filename= file.path(export_fig4, "cfb_tss_cdkn2a.svg"), 
       plot=tss_plot, width = 2.5, height = 5)

tss_plot <- tss_funct("CDKN2B.AS1") +
  labs(x="Treatment", y="Guide z-score")+
  scale_fill_manual(values = c("NT" = "gray90", "CDKN2B.AS1" = "lightskyblue2"), name = "Guide", labels = c("Scrambled control", "CDKN2B-AS1 TSS")) +
  scale_x_discrete(labels = c("Scrambled\ncontrol", "CDKN2B-AS1\nTSS")) +
  stat_compare_means(label = "p.format", method = "t.test", label.x = 1.35, label.y = 1.05) +
  annotate("segment",x = 1.2, xend = 1.8, y = 1, yend = 1)
tss_plot

ggsave(filename= file.path(save_path_tss, "tss_anril.png"), 
       plot=tss_plot, device=png, width = 2.5, height = 5)

ggsave(filename= file.path(export_fig4, "cfb_tss_anril.svg"), 
       plot=tss_plot, width = 2.5, height = 5)
