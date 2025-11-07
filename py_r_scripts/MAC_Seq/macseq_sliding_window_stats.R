#calcualte p-values for MAC-seq sliding windows


library("plotgardener")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(readr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/sliding_window"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/datasets"
input_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/1_platemaps_resources"
save_path_stats = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats"

smc_merged <- read.csv(paste0(save_path_files,"/plate_normalized_guide_means_output.csv"))
#merged_sample_data <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/deseq2/datasets/plate_normalized_guide_means_output.csv")
smc_merged <- smc_merged[,-1]

####make enhancer and SNP plots for use with other graphics####
#enhancer_key <- read.csv(paste0(input_path,"/enhancer_key_controls.csv"))
enhancer_key_only <- read.csv(paste0(input_path,"/enhancer_key.csv"))
#enhancer_key$Start <- as.numeric(gsub(",","",enhancer_key$Start))
#enhancer_key$End <-  as.numeric(gsub(",","",enhancer_key$End))
#enhancer_key$Mid <-  as.numeric(gsub(",","",enhancer_key$Mid))

smc_merged <- read.csv(paste0(save_path_files,"/plate_normalized_guide_means_output.csv"))
#merged_sample_data <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/deseq2/datasets/plate_normalized_guide_means_output.csv")
#smc_merged <- smc_merged[,-1:-4]

#read in permutation data (prm file) and ground truth 500bp binned CRISPRi gene expression (wnd file)
smc_prm_cdkn2b <- read.csv(paste0(save_path_window, "/smc_cdkn2b_permute.csv"))
smc_prm_mtap <- read.csv(paste0(save_path_window, "/smc_mtap_permute.csv"))
smc_prm_dmrta1 <- read.csv(paste0(save_path_window, "/smc_dmrta1_permute.csv"))
smc_prm_cdkn2a <- read.csv(paste0(save_path_window, "/smc_cdkn2a_permute.csv"))
smc_prm_cdkn2b_as1 <- read.csv(paste0(save_path_window, "/smc_cdkn2b_as1_permute.csv"))
smc_prm_mir31hg <- read.csv(paste0(save_path_window, "/smc_mir31hg_permute.csv"))

smc_wnd_cdkn2b <- read.csv(paste0(save_path_window, "/smc_cdkn2b_window.csv"))
smc_wnd_mtap <- read.csv(paste0(save_path_window, "/smc_mtap_window.csv"))
smc_wnd_dmrta1 <- read.csv(paste0(save_path_window, "/smc_dmrta1_window.csv"))
smc_wnd_cdkn2a <- read.csv(paste0(save_path_window, "/smc_cdkn2a_window.csv"))
smc_wnd_cdkn2b_as1 <- read.csv(paste0(save_path_window, "/smc_cdkn2b_as1_window.csv"))
smc_wnd_mir31hg <- read.csv(paste0(save_path_window, "/smc_mir31hg_window.csv"))

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/sliding_window"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/datasets"

cfb_merged <- read.csv(paste0(save_path_files,"/plate_normalized_guide_means_output.csv"))
#merged_sample_data <- read.csv("/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/deseq2/datasets/plate_normalized_guide_means_output.csv")
cfb_merged <- cfb_merged[,-1:-5]

cfb_prm_cdkn2b <- read.csv(paste0(save_path_window, "/cfb_cdkn2b_permute.csv"))
cfb_prm_mtap <- read.csv(paste0(save_path_window, "/cfb_mtap_permute.csv"))
cfb_prm_dmrta1 <- read.csv(paste0(save_path_window, "/cfb_dmrta1_permute.csv"))
cfb_prm_cdkn2a <- read.csv(paste0(save_path_window, "/cfb_cdkn2a_permute.csv"))
cfb_prm_cdkn2b_as1 <- read.csv(paste0(save_path_window, "/cfb_cdkn2b_as1_permute.csv"))
cfb_prm_mir31hg <- read.csv(paste0(save_path_window, "/cfb_mir31hg_permute.csv"))

cfb_wnd_cdkn2b <- read.csv(paste0(save_path_window, "/cfb_cdkn2b_window.csv"))
cfb_wnd_mtap <- read.csv(paste0(save_path_window, "/cfb_mtap_window.csv"))
cfb_wnd_dmrta1 <- read.csv(paste0(save_path_window, "/cfb_dmrta1_window.csv"))
cfb_wnd_cdkn2a <- read.csv(paste0(save_path_window, "/cfb_cdkn2a_window.csv"))
cfb_wnd_cdkn2b_as1 <- read.csv(paste0(save_path_window, "/cfb_cdkn2b_as1_window.csv"))
cfb_wnd_mir31hg <- read.csv(paste0(save_path_window, "/cfb_mir31hg_window.csv"))

####SINGLE REGION CALCULATIONS####
#####CFB MTAP SNP rs1537371#####

#calculate the permutation z-score for an individual window
loc = 22084500
single_perm <- function(loc) {
  x <- cfb_prm_mtap[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(
    prm_mean = mean(perm_val),
    prm_sd = sd(perm_val))
  output <- cfb_wnd_mtap[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (MTAP.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=MTAP.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE),
           gene = "MTAP", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

#calculate MTAP permutation z-score for each position within the locus
cfb_mtap_signif_score_table <- single_perm(22080000)

windows <- seq(22080000,22161500,500)
for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_mtap_signif_score_table <- rbind(cfb_mtap_signif_score_table,z)
}

#single_perm(22084500)
#single_perm(22116000)
#single_perm(22116500)

write.csv(cfb_mtap_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_mtap_window_pvalues.csv")

#hist(x$perm_val,
#     main="Permutation distribution at chr9:22084500 (rs1537371)",
#     xlab="z-score")

####fibroblast CDKN2B permutation z-score and p value calcualtions####

single_perm <- function(loc) {
  x <- cfb_prm_cdkn2b[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(
    prm_mean = mean(perm_val),
    prm_sd = sd(perm_val))
  output <- cfb_wnd_cdkn2b[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2B.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2B.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE),
           gene = "CDKN2B", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

cfb_cdkn2b_signif_score_table <- single_perm(22080000)

windows <- seq(22080000,22161500,500)
for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_cdkn2b_signif_score_table <- rbind(cfb_cdkn2b_signif_score_table,z)
}

write.csv(cfb_cdkn2b_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2b_window_pvalues.csv")

#cdkn2b_signif_score_table <- single_perm(22084500)
#single_perm(22105000)
#single_perm(22105500)

#####SMC MTAP all values#####

single_perm <- function(loc) {
  x <- smc_prm_mtap[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(
    prm_mean = mean(perm_val),
    prm_sd = sd(perm_val))
  output <- smc_wnd_mtap[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (MTAP.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=MTAP.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE),
           gene = "MTAP", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_mtap_signif_score_table <- single_perm(22080000)

windows <- seq(22080000,22161500,500)
for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_mtap_signif_score_table <- rbind(smc_mtap_signif_score_table,z)
}

write.csv(smc_mtap_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_mtap_window_pvalues.csv")


#22115500 #E9
#22116000
#22112000 #E8?



#####SMC CDKN2B all values#####

single_perm <- function(loc) {
  x <- smc_prm_cdkn2b[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(
    prm_mean = mean(perm_val),
    prm_sd = sd(perm_val))
  output <- smc_wnd_cdkn2b[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2B.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2B.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE),
           gene = "CDKN2B", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_cdkn2b_signif_score_table <- single_perm(22080000)

windows <- seq(22080000,22161500,500)
for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_cdkn2b_signif_score_table <- rbind(smc_cdkn2b_signif_score_table,z)
}

write.csv(smc_cdkn2b_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2b_window_pvalues.csv")


#22112000
#22112500

####all other genes: SMC####

single_perm <- function(loc) {
  x <- smc_prm_cdkn2a[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- smc_wnd_cdkn2a[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2A.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2A.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "CDKN2A", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_cdkn2a_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_cdkn2a_signif_score_table <- rbind(smc_cdkn2a_signif_score_table,z)
}

write.csv(smc_cdkn2a_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2a_window_pvalues.csv")

single_perm <- function(loc) {
  x <- smc_prm_cdkn2b_as1[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- smc_wnd_cdkn2b_as1[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2B.AS1.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2B.AS1.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "CDKN2B.AS1", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_cdkn2b_as1_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_cdkn2b_as1_signif_score_table <- rbind(smc_cdkn2b_as1_signif_score_table,z)
}

write.csv(smc_cdkn2b_as1_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2b_as1_window_pvalues.csv")

single_perm <- function(loc) {
  x <- smc_prm_dmrta1[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- smc_wnd_dmrta1[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (DMRTA1.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=DMRTA1.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "DMRTA1", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_dmrta1_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_dmrta1_signif_score_table <- rbind(smc_dmrta1_signif_score_table,z)
}

write.csv(smc_dmrta1_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_dmrta1_window_pvalues.csv")

single_perm <- function(loc) {
  x <- smc_prm_mir31hg[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- smc_wnd_mir31hg[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (MIR31HG.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=MIR31HG.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "MIR31HG", celltype = "SMC")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

smc_mir31hg_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  smc_mir31hg_signif_score_table <- rbind(smc_mir31hg_signif_score_table,z)
}

write.csv(smc_mir31hg_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_mir31hg_window_pvalues.csv")

####all other genes: CFB####

single_perm <- function(loc) {
  x <- cfb_prm_cdkn2a[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- cfb_wnd_cdkn2a[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2A.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2A.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "CDKN2A", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

cfb_cdkn2a_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_cdkn2a_signif_score_table <- rbind(cfb_cdkn2a_signif_score_table,z)
}

write.csv(cfb_cdkn2a_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2a_window_pvalues.csv")

single_perm <- function(loc) {
  x <- cfb_prm_cdkn2b_as1[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- cfb_wnd_cdkn2b_as1[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (CDKN2B.AS1.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=CDKN2B.AS1.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "CDKN2B.AS1", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

cfb_cdkn2b_as1_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_cdkn2b_as1_signif_score_table <- rbind(cfb_cdkn2b_as1_signif_score_table,z)
}

write.csv(cfb_cdkn2b_as1_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2b_as1_window_pvalues.csv")

single_perm <- function(loc) {
  x <- cfb_prm_dmrta1[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- cfb_wnd_dmrta1[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (DMRTA1.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=DMRTA1.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "DMRTA1", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

cfb_dmrta1_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_dmrta1_signif_score_table <- rbind(cfb_dmrta1_signif_score_table,z)
}

write.csv(cfb_dmrta1_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_dmrta1_window_pvalues.csv")

single_perm <- function(loc) {
  x <- cfb_prm_mir31hg[-1] %>% dplyr::filter(coords == loc)
  y <- x %>% summarize(prm_mean = mean(perm_val),prm_sd = sd(perm_val))
  output <- cfb_wnd_mir31hg[-1] %>% dplyr::filter(coords == loc) %>% cbind(y) %>%
    mutate(zscore = (MIR31HG.mean-prm_mean)/prm_sd,
           pval = 1*pnorm(q=MIR31HG.mean, mean = prm_mean, sd = prm_sd, lower.tail=TRUE), 
           gene = "MIR31HG", celltype = "CFB")
  colnames(output)[1] <- "norm_expression"
  return(output)
}

cfb_mir31hg_signif_score_table <- single_perm(22080000)

for (i in 2:164) {
  z <- single_perm(windows[i])
  cfb_mir31hg_signif_score_table <- rbind(cfb_mir31hg_signif_score_table,z)
}

write.csv(cfb_mir31hg_signif_score_table, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_mir31hg_window_pvalues.csv")


####save permuted values####

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/sliding_window"

smc_wnd_cdkn2b <- cbind(smc_wnd_cdkn2b[-1],smc_cdkn2b_signif_score_table[4:8])
smc_wnd_mtap <- cbind(smc_wnd_mtap[-1],smc_mtap_signif_score_table[4:8])
smc_wnd_dmrta1 <- cbind(smc_wnd_dmrta1[-1],smc_dmrta1_signif_score_table[4:8])
smc_wnd_cdkn2a <- cbind(smc_wnd_cdkn2a[-1],smc_cdkn2a_signif_score_table[4:8])
smc_wnd_cdkn2b_as1 <- cbind(smc_wnd_cdkn2b_as1[-1],smc_cdkn2b_as1_signif_score_table[4:8])
smc_wnd_mir31hg <- cbind(smc_wnd_mir31hg[-1],smc_mir31hg_signif_score_table[4:8])

smc_wnd_cdkn2b <- smc_wnd_cdkn2b %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
smc_wnd_mtap <- smc_wnd_mtap %>% dplyr::select(coords, everything()) %>%mutate(n_guides = floor(n_guides/2))
smc_wnd_dmrta1 <- smc_wnd_dmrta1 %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
smc_wnd_cdkn2a <- smc_wnd_cdkn2a %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
smc_wnd_cdkn2b_as1 <- smc_wnd_cdkn2b_as1 %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
smc_wnd_mir31hg <- smc_wnd_mir31hg %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))

colnames(smc_wnd_cdkn2b)[2] <- "norm_expression"
colnames(smc_wnd_mtap)[2] <- "norm_expression"
colnames(smc_wnd_dmrta1)[2] <- "norm_expression"
colnames(smc_wnd_cdkn2a)[2] <- "norm_expression"
colnames(smc_wnd_cdkn2b_as1)[2] <- "norm_expression"
colnames(smc_wnd_mir31hg)[2] <- "norm_expression"

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/sliding_window"

cfb_wnd_cdkn2b <- cbind(cfb_wnd_cdkn2b[-1],cfb_cdkn2b_signif_score_table[4:8])
cfb_wnd_mtap <- cbind(cfb_wnd_mtap[-1],cfb_mtap_signif_score_table[4:8])
cfb_wnd_dmrta1 <- cbind(cfb_wnd_dmrta1[-1],cfb_dmrta1_signif_score_table[4:8])
cfb_wnd_cdkn2a <- cbind(cfb_wnd_cdkn2a[-1],cfb_cdkn2a_signif_score_table[4:8])
cfb_wnd_cdkn2b_as1 <- cbind(cfb_wnd_cdkn2b_as1[-1],cfb_cdkn2b_as1_signif_score_table[4:8])
cfb_wnd_mir31hg <- cbind(cfb_wnd_mir31hg[-1],cfb_mir31hg_signif_score_table[4:8])

cfb_wnd_cdkn2b <- cfb_wnd_cdkn2b %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
cfb_wnd_mtap <- cfb_wnd_mtap %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
cfb_wnd_dmrta1 <- cfb_wnd_dmrta1 %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
cfb_wnd_cdkn2a <- cfb_wnd_cdkn2a %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2)) 
cfb_wnd_cdkn2b_as1 <- cfb_wnd_cdkn2b_as1 %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))
cfb_wnd_mir31hg <- cfb_wnd_mir31hg %>% dplyr::select(coords, everything()) %>% mutate(n_guides = floor(n_guides/2))

colnames(cfb_wnd_cdkn2b)[2] <- "norm_expression"
colnames(cfb_wnd_mtap)[2] <- "norm_expression"
colnames(cfb_wnd_dmrta1)[2] <- "norm_expression"
colnames(cfb_wnd_cdkn2a)[2] <- "norm_expression"
colnames(cfb_wnd_cdkn2b_as1)[2] <- "norm_expression"
colnames(cfb_wnd_mir31hg)[2] <- "norm_expression"

####save data####

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats"

write.csv(smc_wnd_cdkn2b, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2b_window_pvalues.csv")
write.csv(smc_wnd_mtap, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_mtap_window_pvalues.csv")
write.csv(smc_wnd_dmrta1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_dmrta1_window_pvalues.csv")
write.csv(smc_wnd_cdkn2b_as1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2b_as1_window_pvalues.csv")
write.csv(smc_wnd_cdkn2a, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_cdkn2a_window_pvalues.csv")
write.csv(smc_wnd_mir31hg, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/perm_stats/smc_mir31hg_window_pvalues.csv")

save_path_window = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/sliding_window"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/datasets"

write.csv(cfb_wnd_cdkn2b, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2b_window_pvalues.csv")
write.csv(cfb_wnd_mtap, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_mtap_window_pvalues.csv")
write.csv(cfb_wnd_dmrta1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_dmrta1_window_pvalues.csv")
write.csv(cfb_wnd_cdkn2b_as1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2b_as1_window_pvalues.csv")
write.csv(cfb_wnd_cdkn2a, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_cdkn2a_window_pvalues.csv")
write.csv(cfb_wnd_mir31hg, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/cfb/zscore_normalized/perm_stats/cfb_mir31hg_window_pvalues.csv")


####test FDR
test <- rbind(cfb_wnd_cdkn2b,cfb_wnd_mtap,cfb_wnd_dmrta1,cfb_wnd_cdkn2b_as1,cfb_wnd_cdkn2a,cfb_wnd_mir31hg)

hist(test$pval,breaks = 20,
     main="P-value distribution for all genes",
     xlab="z-score")

hist(cfb_wnd_mtap$pval,breaks = 20,
     main="P-value distribution for MTAP",
     xlab="z-score")
