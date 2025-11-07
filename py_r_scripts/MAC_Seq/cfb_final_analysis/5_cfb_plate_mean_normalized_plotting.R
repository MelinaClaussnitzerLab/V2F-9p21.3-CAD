#

library(dplyr)
library(ggplot2)
library(reghelper)
library(patchwork)
library(ggrepel)
library(svglite)
library(tidyverse)
library(ggpubr)

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

#create simple enhancer graphic for combination with other plots
enhancer_locations <- ggplot(data=enhancer_key) +
  xlim(22064500,22135000) + #xlim(22064500,22163000) + ylim (-1,1) +
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
#enhancer_plot

#create simple enhancer graphic for combination with other plots
snp_locations <- ggplot(data=snp_key) +
  geom_rect(data = snp_key, inherit.aes = F, aes(
    xmin = loc-10, xmax = loc+10, ymin = -0.05,ymax = 0.05, fill=is.na(constraint), color=is.na(constraint))) +
  xlim(22064500,22135000) + #xlim(22064500,22163000) + #ylim (-1,1) +
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
#snp_plot

###TESTING####

#a=ggplot object
#b=coefs
#sd(cfb_merged$CDKN2A)
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
               aes(x=GuidePosition,y=h[[g]],color="red")) +
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

####use sliding window approach to look at average CRISPRi signal####

windows <- seq(22080000,22161750,500)

# z <- list()
# q <- 1
# for (i in windows) {
#   x <- cfb_merged %>% filter((i-500)<GuidePosition & (i+500)>GuidePosition)
#   y <- mean(x$CDKN2B.mean)
#   z[q] <- y
#   q <- q+1
# }

sliding_window <- function(a,b) {
  z <- list()
  windows <- seq(22080000,22161750,500)
  q <- 1
  for (i in windows) {
    x <- a %>% filter((i-500)<GuidePosition & (i+500)>GuidePosition)
    y <- mean(x[[b]])
    z[q] <- y
    q <- q+1
  }
  c <- t(data.frame(z))
  c <- data.frame(c)
  rownames(c) <- windows
  colnames(c) <- b
  return(c)
}

sliding_window_n <- function(a,b) {
  z <- list()
  s <- list()
  windows <- seq(22080000,22161750,500)
  q <- 1
  for (i in windows) {
    u <- a %>% filter((i-500)<GuidePosition & (i+500)>GuidePosition)
    r <- mean(u[[b]])
    s[q] <- nrow(u)
    z[q] <- r
    q <- q+1
  }
  c <- t(data.frame(z))
  rownames(c) <- windows
  colnames(c) <- b
  d <- t(data.frame(s))
  rownames(d) <- windows
  colnames(d) <- "n_guides"
  e <- data.frame(cbind(c,d))
  e <- e %>% mutate(coords = as.numeric(row.names(e)))
  return(e)
}

test <- c()
test <-  sliding_window_n(cfb_merged, "CDKN2A.mean")

cdkn2b_window <-  sliding_window_n(cfb_merged, "CDKN2B.mean")

cdkn2a_window <-  sliding_window_n(cfb_merged, "CDKN2A.mean")

anril_window <-  sliding_window_n(cfb_merged, "CDKN2B.AS1.mean")

dmrta1_window <-  sliding_window_n(cfb_merged, "DMRTA1.mean")

mtap_window <-  sliding_window_n(cfb_merged, "MTAP.mean")

mir31hg_window <-  sliding_window_n(cfb_merged, "MIR31HG.mean")

cdkn2b_plot <- nnfplot(cdkn2b_plot,cfb_merged, "CDKN2B.mean")
y <- cdkn2b_plot + ggtitle("CDKN2B") + loc_lims +
  geom_line(data=cdkn2b_window, aes(x=coords,y=CDKN2B.mean, color="red"))
y / snp_locations / enhancer_locations

ggsave(filename= file.path(save_path_window, "cdkn2b_sliding_window.png"), 
       plot=y, device=png, width = 10, height = 5)

anril_plot <- nnfplot(anril_plot,cfb_merged, "CDKN2B.AS1.mean")
y <- anril_plot + ggtitle("ANRIL") + loc_lims +
  geom_line(data=anril_window, aes(x=coords,y=CDKN2B.AS1.mean, color="red"))
y
ggsave(filename= file.path(save_path_window, "anril_sliding_window.png"), 
       plot=y, device=png, width = 10, height = 5)

cdkn2a_plot <- nnfplot(cdkn2a_plot,cfb_merged, "CDKN2A.mean")
y <- cdkn2a_plot + ggtitle("CDKN2A") + loc_lims + 
  geom_line(data=cdkn2a_window, aes(x=coords,y=CDKN2A.mean, color="red"))
y
ggsave(filename= file.path(save_path_window, "cdkn2a_sliding_window.png"), 
       plot=y, device=png, width = 10, height = 5)

dmrta1_plot <- nnfplot(dmrta1_plot,cfb_merged, "DMRTA1.mean") #+ xlim(22064500,22135000) #xlim(22064500,22163000)
w <- dmrta1_plot +  ggtitle("DMRTA1") + loc_lims +
  geom_line(data=dmrta1_window, aes(x=coords,y=DMRTA1.mean, color="red"))
w
ggsave(filename= file.path(save_path_window, "dmrta1_sliding_window.png"), 
       plot=w, device=png, width = 10, height = 5)

mtap_plot <- nnfplot(mtap_plot,cfb_merged, "MTAP.mean") #+ xlim(22064500,22135000) #xlim(22064500,22163000)
v <- mtap_plot + ggtitle("MTAP") + loc_lims +
  geom_line(data=mtap_window, aes(x=as.numeric(rownames(mtap_window)),y=MTAP.mean, color="red"))
v
ggsave(filename= file.path(save_path_window, "mtap_sliding_window.png"), 
       plot=v, device=png, width = 10, height = 5)
y <- v/enhancer_locations/snp_locations + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
y

mir31hg_plot <- nnfplot(mtap_plot,cfb_merged, "MIR31HG.mean") #+ xlim(22064500,22135000) #xlim(22064500,22163000)
v <- mtap_plot + ggtitle("MIR31HG") + loc_lims +
  geom_line(data=mir31hg_window, aes(x=as.numeric(rownames(mir31hg_window)),y=MIR31HG.mean, color="red"))
v
ggsave(filename= file.path(save_path_window, "mir31hg_sliding_window.png"), 
       plot=v, device=png, width = 10, height = 5)
y <- v/enhancer_locations/snp_locations + plot_layout(ncol = 1, nrow=3, heights = c(1, 0.3,0.3))
y

####run permutation testing####

simple_plot <- function(a,b,g) {
  a <- ggplot(b, aes(x=GuidePosition, y=b[[g]])) +
    geom_point(color = "darkgray") +
    #geom_point(data=h, inherit.aes = F, 
    #           aes(x=GuidePosition,y=h[[g]],color="red")) +
    geom_hline(yintercept= 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept= -1, linetype = "dashed", color = "black") +
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
}

# cdkn2b_permute <- cdkn2b_window
# cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
# n_buckets <- cdkn2b_permute$n_guides
# 
# for (q in 1:500) {
#   temp <- data.frame(rep(0,164))
#   colnames(temp) <- "perm_val"
#   temp <- temp %>% mutate(perm_rep = paste0("V",q))
#   #set.seed(q)
#   for (i in 1:164) {
#     temp[i,1] <- mean(sample(cfb_sample$CDKN2B.mean, n_buckets[[i]]))
#   }
#   #temp <- as.data.frame(t(temp))
#   #test1 <- cbind(test1,temp)
#   temp_var <- cbind(cdkn2b_window, temp)
#   if (q==1){
#     cdkn2b_permute <- temp_var
#   } else if (q!=1){
#     cdkn2b_permute <- rbind(cdkn2b_permute, temp_var)
#   }
# }

cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
cfb_sample_unique <- cfb_sample[!duplicated(cfb_sample$GuidePosition),]
#test_window <-  sliding_window_n(cfb_sample_unique, "CDKN2B.mean")
for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$CDKN2B.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    cdkn2b_permute <- temp
  } else if (q!=1){
    cdkn2b_permute <- rbind(cdkn2b_permute, temp)
  }
}

cdkn2b_plot <- simple_plot(cdkn2b_plot,cfb_merged, "CDKN2B.mean")
y <- cdkn2b_plot + ggtitle("CDKN2B") + loc_lims + guides(color="none") +
  geom_line(data=cdkn2b_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=cdkn2b_window, aes(x=as.numeric(rownames(cdkn2b_window)),y=CDKN2B.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "cdkn2b_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)


#permute MTAP values
#mtap_permute <- mtap_window
#cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
#n_buckets <- mtap_permute$n_guides

#test_window <-  sliding_window_n(cfb_sample_unique, "MTAP.mean")
for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$MTAP.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    mtap_permute <- temp
  } else if (q!=1){
    mtap_permute <- rbind(mtap_permute, temp)
  }
}

mtap_plot <- simple_plot(mtap_plot,cfb_merged, "MTAP.mean")
y <- mtap_plot + ggtitle("MTAP") + loc_lims + guides(color="none") +
  geom_line(data=mtap_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=mtap_window, aes(x=coords,y=MTAP.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "mtap_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)

#permute CDKN2A values
#cdkn2a_permute <- cdkn2a_window
#cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
#n_buckets <- cdkn2a_permute$n_guides

#test_window <-  sliding_window_n(cfb_sample_unique, "MTAP.mean")
for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$CDKN2A.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    cdkn2a_permute <- temp
  } else if (q!=1){
    cdkn2a_permute <- rbind(cdkn2a_permute, temp)
  }
}

cdkn2a_plot <- simple_plot(cdkn2a_plot,cfb_merged, "CDKN2A.mean")
y <- mtap_plot + ggtitle("MTAP") + loc_lims + guides(color="none") +
  geom_line(data=cdkn2a_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=cdkn2a_window, aes(x=coords,y=CDKN2A.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "cdkn2a_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)

#permute ANRIL values
#anril_permute <- anril_window
#cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
#n_buckets <- anril_permute$n_guides

for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$CDKN2B.AS1.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    anril_permute <- temp
  } else if (q!=1){
    anril_permute <- rbind(anril_permute, temp)
  }
}

anril_plot <- simple_plot(anril_plot,cfb_merged, "CDKN2B.AS1.mean")
y <- mtap_plot + ggtitle("CDKN2B-AS1") + loc_lims + guides(color="none") +
  geom_line(data=anril_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=anril_window, aes(x=coords,y=CDKN2B.AS1.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "anril_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)


#permute DMRTA1 values
#dmrta1_permute <- dmrta1_window
#cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
#n_buckets <- dmrta1_permute$n_guides

for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$DMRTA1.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    dmrta1_permute <- temp
  } else if (q!=1){
    dmrta1_permute <- rbind(dmrta1_permute, temp)
  }
}

dmrta1_plot <- simple_plot(dmrta1_plot,cfb_merged, "DMRTA1.mean")
y <- dmrta1_plot + ggtitle("DMRTA1") + loc_lims + guides(color="none") +
  geom_line(data=dmrta1_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=dmrta1_window, aes(x=coords,y=DMRTA1.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "dmrta1_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)


y / snp_locations / enhancer_locations

#permute MIR31HG values
mir31hg_permute <- mir31hg_window
cfb_sample <- cfb_merged %>% filter(guide_ctrl == "Screen")
n_buckets <- mir31hg_permute$n_guides

for (q in 1:500) {
  x <- cfb_sample_unique
  x$perm_val <- sample(cfb_sample_unique$MIR31HG.mean, 508)
  temp <-  sliding_window_n(x, "perm_val")
  temp <- temp %>% mutate(perm_rep = paste0("V",q))
  #set.seed(q)
  #temp <- as.data.frame(t(temp))
  #test1 <- cbind(test1,temp)
  if (q==1){
    mir31hg_permute <- temp
  } else if (q!=1){
    mir31hg_permute <- rbind(mir31hg_permute, temp)
  }
}

mir31hg_plot <- simple_plot(mir31hg_plot,cfb_merged, "MIR31HG.mean")
y <- mir31hg_plot + ggtitle("MIR31HG") + loc_lims + guides(color="none") +
  geom_line(data=mir31hg_permute, aes(x=coords,y=perm_val, color=perm_rep), alpha = 0.1) +
  geom_line(data=mir31hg_window, aes(x=coords,y=MIR31HG.mean), color="black")
y

ggsave(filename= file.path(save_path_window, "mir31hg_permutation.png"), 
       plot=y, device=png, width = 8, height = 5)


####save all values####
write.csv(cdkn2b_permute, paste0(save_path_window, "/cfb_cdkn2b_permute.csv"))
write.csv(mtap_permute, paste0(save_path_window, "/cfb_mtap_permute.csv"))
write.csv(dmrta1_permute, paste0(save_path_window, "/cfb_dmrta1_permute.csv"))
write.csv(cdkn2a_permute, paste0(save_path_window, "/cfb_cdkn2a_permute.csv"))
write.csv(anril_permute, paste0(save_path_window, "/cfb_cdkn2b_as1_permute.csv"))
write.csv(mir31hg_permute, paste0(save_path_window, "/cfb_mir31hg_permute.csv"))

write.csv(cdkn2b_window, paste0(save_path_window, "/cfb_cdkn2b_window.csv"))
write.csv(mtap_window, paste0(save_path_window, "/cfb_mtap_window.csv"))
write.csv(dmrta1_window, paste0(save_path_window, "/cfb_dmrta1_window.csv"))
write.csv(cdkn2a_window, paste0(save_path_window, "/cfb_cdkn2a_window.csv"))
write.csv(anril_window, paste0(save_path_window, "/cfb_cdkn2b_as1_window.csv"))
write.csv(mir31hg_window, paste0(save_path_window, "/cfb_mir31hg_window.csv"))
####plot sliding window correlation by replicate####

r2_funct <- function(z) {
  model <-summary(lm(rep1~rep2,data=z))
  r2 = round(model$r.squared,3)
  substitute(italic(r)^2~"="~r2)
}

simple_theme <- theme_classic() +
  theme(text = element_text(size = 15),
        title = element_text(size = 20),
        legend.text = element_text(size=15),
        axis.title.x=element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15))

#generate replicates

cfb_r1 <- cfb_merged %>% filter(rep_id %% 2 == 0)
cfb_r2 <- cfb_merged %>% filter(rep_id %% 2 == 1)


window_reps <- sliding_window(cfb_r1, "CDKN2B")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "CDKN2B")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

cdkn2b_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey")+
  annotate(geom="text", x=-1, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "CDKN2B sliding window agreement")
ggsave(filename= file.path(save_path_window, "cdkn2b_window_correlation.png"), 
       plot=cdkn2b_plot_agreement, device=png, width = 5, height = 5)

window_reps <- sliding_window(cfb_r1, "CDKN2A")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "CDKN2A")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

cdkn2a_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey")+
  annotate(geom="text", x=-1.5, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "CDKN2A sliding window agreement")
ggsave(filename= file.path(save_path_window, "cdkn2a_window_correlation.png"), 
       plot=cdkn2a_plot_agreement, device=png, width = 5, height = 5)

window_reps <- sliding_window(cfb_r1, "MTAP")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "MTAP")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

mtap_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey") +
  annotate(geom="text", x=-1.5, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "MTAP sliding window agreement")
ggsave(filename= file.path(save_path_window, "mtap_window_correlation.png"), 
       plot=mtap_plot_agreement, device=png, width = 5, height = 5)

window_reps <- sliding_window(cfb_r1, "CDKN2B.AS1")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "CDKN2B.AS1")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

anril_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey") +
  annotate(geom="text", x=-1, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "CDKN2B-AS1 sliding window agreement")
ggsave(filename= file.path(save_path_window, "anril_window_correlation.png"), 
       plot=anril_plot_agreement, device=png, width = 5, height = 5)

window_reps <- sliding_window(cfb_r1, "DMRTA1")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "DMRTA1")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

dmrta1_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey") +
  annotate(geom="text", x=-0.5, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "DMRTA1 sliding window agreement")
ggsave(filename= file.path(save_path_window, "dmrta1_window_correlation.png"), 
       plot=dmrta1_plot_agreement, device=png, width = 5, height = 5)

mir31hg_reps <- sliding_window(cfb_r1, "MIR31HG")
rownames(window_reps) <-  windows
temp <- sliding_window(cfb_r2, "MIR31HG")
rownames(temp) <- windows
window_reps <- data.frame(cbind(window_reps,temp))
colnames(window_reps) <- c("rep1","rep2")

mir31hg_plot_agreement <- ggplot(window_reps, aes(x=rep1,y=rep2)) +
  geom_point(alpha=0.6) + geom_smooth(method = lm, color="grey") +
  annotate(geom="text", x=-0.5, y=1, label = as.call(r2_funct(window_reps))) +
  simple_theme + labs(x = "Rep 1 sliding window value", y = "Rep 2 sliding window value", title = "DMRTA1 sliding window agreement")
ggsave(filename= file.path(save_path_window, "mir31hg_window_correlation.png"), 
       plot=mir31hg_plot_agreement, device=png, width = 5, height = 5)

#model_mtap <-summary(lm(rep1~rep2,data=mtap_window_reps))
