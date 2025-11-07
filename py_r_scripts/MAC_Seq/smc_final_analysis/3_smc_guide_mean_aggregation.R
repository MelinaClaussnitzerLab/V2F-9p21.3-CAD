#Script for grouping MAC-Seq guides by mean for further plotting

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
guide_input_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/input"
save_path_figs = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/locus_plots"
save_path_files = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/datasets"
save_path_batch = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/whole_screen_output/smc/zscore_normalized/batch_effect"

smc_merged <- read.csv(paste0(guide_input_path,"/plate_normalized_cpm_reads_all_input.csv"))

merged_means_calc <- smc_merged %>% 
  select(c(guide_name,guide_id,CDKN2B.AS1,MTAP,CDKN2A,CDKN2B,DMRTA1, MIR31HG)) %>% 
  group_by(guide_id) %>% 
  mutate(
    CDKN2B.AS1.mean = mean(CDKN2B.AS1),
    MTAP.mean = mean(MTAP),
    CDKN2A.mean = mean(CDKN2A),
    CDKN2B.mean = mean(CDKN2B),
    DMRTA1.mean = mean(DMRTA1),
    MIR31HG.mean = mean(MIR31HG)
  ) %>%
  select(c(guide_name,guide_id,CDKN2B.AS1.mean,MTAP.mean,CDKN2A.mean,CDKN2B.mean,DMRTA1.mean, MIR31HG.mean))

smc_merged <- merge(smc_merged,merged_means_calc,by=c("guide_name","guide_id"))

smc_merged <- smc_merged %>% mutate(
  guide_ctrl = ifelse(
    guide == "MTAP" | guide == "CDKN2A" | guide == "CDKN2B" | guide == "CDKN2B-AS1" |
      guide == "DMRTA1" | guide == "CDKN2A.DT", "TSS",
    guide_ctrl
  )
)

####important enhancer database and annotate nearby enhancers####
enhancer_key <- read.csv(paste0(input_path,"/enhancer_key_controls.csv"))
#enhancer_key_only <- read.csv("raw_data/enhancer_key.csv")
enhancer_key_only <- read.csv(paste0(input_path,"/enhancer_key.csv"))
enhancer_key$Start <- as.numeric(gsub(",","",enhancer_key$Start))
enhancer_key$End <-  as.numeric(gsub(",","",enhancer_key$End))
enhancer_key$Mid <-  as.numeric(gsub(",","",enhancer_key$Mid))

#define graphed position for negative controls (2 kinds), and assign these locations in the data
####assign plot positions for negstive controls####
negguidepos = 22065000
locguidepos = 22070000


#label all controls with the artificial "genomic position" they will be graphed on
smc_merged <- smc_merged %>%
  mutate(GuidePositionOrig = GuidePosition,
         GuidePosition = ifelse(guide == "NEG_CONTROL", negguidepos,
                                ifelse(guide == "Neg-9p21 regions", locguidepos,
                                       GuidePosition)),
         mtper = mtper*0.01,
         riboper = riboper*0.01) %>%
  relocate(GuidePositionOrig, .after = GuidePosition)

enhancer_samples <- c("")
for (i in 1:nrow(enhancer_key)) {
  z <- enhancer_key[i,2]
  x <- enhancer_key[i,3]
  y <- smc_merged %>% filter(z<GuidePosition & x>GuidePosition)
  y <- y %>% mutate(
    near = enhancer_key[i,1])
  enhancer_samples <- rbind(enhancer_samples,y)
}
enhancer_samples <- enhancer_samples[-1,]
enhancer_samples <- enhancer_samples %>% select(guide_name,near)

smc_merged <- merge(smc_merged,enhancer_samples, by="guide_name", all.x=T)

smc_merged <- smc_merged %>% mutate(
  near = ifelse(
    guide == "MTAP", "MTAP", ifelse(guide == "CDKN2A", "CDKN2A", ifelse(guide == "CDKN2B", "CDKN2B",
                                                                        ifelse(guide == "CDKN2B-AS1", "CDKN2B.AS1", ifelse(guide == "DMRTA1", "DMRTA1", 
                                                                                                                           ifelse(guide == "CDKN2A.DT", "CDKN2A.DT", ifelse(is.na(near)==T, "None",
                                                                                                                                                                            near))))))
  ))

#sort sample data by metadata for a consistent order
smc_merged <- smc_merged %>% arrange(plate,num)
#rownames(smc_merged) <-  smc_merged$guide_name

smc_merged_dedup <- smc_merged[!duplicated(smc_merged$guide_id), ]
smc_merged_dedup <- smc_merged_dedup[-23:-28]
write.csv(smc_merged_dedup, file = paste0(save_path_files, "/plate_normalized_guide_means_deduplicated.csv"))

write.csv(smc_merged, file = paste0(save_path_files, "/plate_normalized_guide_means_output.csv"))

