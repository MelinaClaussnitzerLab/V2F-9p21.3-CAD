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








#-----------effect of treat
treatmentVSControl <- function(data){
  
  results <- list()
  idx <- 1
  
  meta_cols <- grep("^Metadata_", names(data), value = TRUE)
  
  for (plate in unique(data$Metadata_Plate)) {
    plate_df <- data %>% filter(Metadata_Plate == plate)
    print(plate)
    
    for (stim in unique(plate_df$Metadata_stimulation)) {
      print(stim)
      stim_df <- plate_df %>% filter(Metadata_stimulation == stim)
      
      ref_df <- stim_df %>% filter(Metadata_treatment == "nontargeting control")
      
      feature_cols <- setdiff(names(stim_df), meta_cols)
      
      other_treats <- setdiff(unique(stim_df$Metadata_treatment),
                              "nontargeting control")
      for (treatment in other_treats) {
        treat_df <- stim_df %>% filter(Metadata_treatment == treatment)
        
        for (feat in feature_cols) {
          ref_vals <- as.numeric(ref_df[[feat]])
          tr_vals  <- as.numeric(treat_df[[feat]])
          
          ref_mean <- mean(ref_vals, na.rm = TRUE)
          tr_mean  <- mean(tr_vals,  na.rm = TRUE)
          ref_n    <- sum(!is.na(ref_vals))
          tr_n     <- sum(!is.na(tr_vals))
          direction <- ifelse(tr_mean >= ref_mean, "inc", "dec")
          
          tryCatch({
            wtest <- suppressWarnings(
              wilcox.test(ref_vals, tr_vals, conf.int = TRUE)
            )
            es <- wilcox_effsize(
              data.frame(
                value    = c(ref_vals, tr_vals),
                category = rep(c("reference","treatment"),
                               c(length(ref_vals), length(tr_vals)))
              ),
              value ~ category
            )$effsize
            
            results[[idx]] <- data.frame(
              reference       = "nontargeting control",
              treatment       = treatment,
              stimulation     = stim,
              pvalue          = wtest$p.value,
              Metadata_Plate  = plate,
              featureName     = feat,
              effectSize      = es,
              direction       = direction,
              referenceMean   = ref_mean,
              treatmentMean   = tr_mean,
              referenceN      = ref_n,
              treatmentN      = tr_n,
              stringsAsFactors = FALSE
            )
            idx <- idx + 1
          }, error = function(e){
            message(sprintf(
              "Wilcox test failed for feature '%s' on plate '%s', stimulation '%s', treatment '%s' vs reference: %s",
              feat, plate, stim, treatment, e$message
            ))
          })
        }
      }
    }
  }
  
  if (length(results) > 0) {
    tdf <- do.call(rbind, results)
  } else {
    tdf <- data.frame(
      reference      = character(),
      treatment      = character(),
      stimulation    = character(),
      pvalue         = numeric(),
      Metadata_Plate = character(),
      featureName    = character(),
      effectSize     = numeric(),
      direction      = character(),
      referenceMean  = numeric(),
      treatmentMean  = numeric(),
      referenceN     = integer(),
      treatmentN     = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  return(tdf)
}



smc <- read.csv('/Users/cho/Documents/cellprofiler/downstream/9p21/072225/int_smc.csv')
head(smc)


tdf_smc <- treatmentVSControl(smc)
head(tdf_smc)



tdf <- tdf_smc %>%
  group_by(treatment, stimulation) %>%
  mutate(qvalue = qvalue(pvalue)$qvalues) %>%
  ungroup()


tdf <- tdf %>%
  arrange(Metadata_Plate,stimulation,treatment,qvalue)

write.csv(tdf, "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/smc_effectOfKD.csv", row.names = FALSE)













# basalVStgfb <- function(data){
#   results <- list()
#   idx <- 1
#   
#   meta_cols <- grep("^Metadata_", names(data), value = TRUE)
#   
#   for(plate in unique(data$Metadata_Plate)){
#     plate_df <- data %>% filter(Metadata_Plate == plate)
#     print(plate)
#     
#     for(treatment in unique(plate_df$Metadata_treatment)){
#       print(treatment)
#       tr_df <- plate_df %>% filter(Metadata_treatment == treatment)
#       
#       basal_df <- tr_df %>% filter(Metadata_stimulation == "basal")
#       tgfb_df  <- tr_df %>% filter(Metadata_stimulation == "tgfb")
#       
#       feature_cols <- setdiff(names(tr_df), meta_cols)
#       
#       for(feat in feature_cols){
#         b_vals <- as.numeric(basal_df[[feat]])
#         t_vals <- as.numeric(tgfb_df[[feat]])
#         
#         b_mean <- mean(b_vals, na.rm = TRUE)
#         t_mean <- mean(t_vals, na.rm = TRUE)
#         b_n    <- sum(!is.na(b_vals))
#         t_n    <- sum(!is.na(t_vals))
#         direction <- ifelse(t_mean >= b_mean, "inc", "dec")
#         
#         tryCatch({
#           wtest <- suppressWarnings(wilcox.test(b_vals, t_vals, conf.int = T))
#           es <- wilcox_effsize(
#             data.frame(
#               value    = c(b_vals, t_vals),
#               category = rep(c("basal","tgfb"), c(length(b_vals), length(t_vals)))
#             ),
#             value ~ category
#           )$effsize
#           
#           results[[idx]] <- data.frame(
#             basal        = treatment,
#             tgfb         = treatment,
#             pvalue       = wtest$p.value,
#             Metadata_Plate = plate,
#             featureName  = feat,
#             effectSize   = es,
#             direction    = direction,
#             basalMean    = b_mean,
#             tgfbMean     = t_mean,
#             basalN       = b_n,
#             tgfbN        = t_n,
#             stringsAsFactors = FALSE
#           )
#           idx <- idx + 1
#         }, error = function(e){
#           message(
#             sprintf(
#               "Wilcox test failed for feature '%s' on plate '%s', treatment '%s': %s",
#               feat, plate, treatment, e$message
#             )
#           )
#         })
#       }
#     }
#   }
#   
#   if(length(results) > 0){
#     tdf <- do.call(rbind, results)
#   } else {
#     tdf <- data.frame(
#       basal=character(), tgfb=character(), pvalue=numeric(),
#       Metadata_Plate=character(), featureName=character(),
#       effectSize=numeric(), direction=character(),
#       basalMean=numeric(), tgfbMean=numeric(),
#       basalN=integer(), tgfbN=integer(),
#       stringsAsFactors = FALSE
#     )
#   }
#   
#   return(tdf)
# }
# 
# 
# 
# 
# smc <- read.csv('/Users/cho/Documents/cellprofiler/downstream/9p21/072225/int_smc.csv')
# head(smc)
# 
# 
# tdf_smc <- basalVStgfb(smc)
# head(tdf_smc)
# 
# 
# 
# tdf <- tdf_smc %>%
#   group_by(basal) %>%
#   mutate(qvalue = qvalue(pvalue)$qvalues) %>%
#   ungroup()
# 
# 
# tdf <- tdf %>%
#   arrange(Metadata_Plate,basal, qvalue)
# 
# write.csv(tdf, "/Users/cho/Documents/cellprofiler/downstream/9p21/072225/smc_effectOfStim.csv", row.names = FALSE)
# 










