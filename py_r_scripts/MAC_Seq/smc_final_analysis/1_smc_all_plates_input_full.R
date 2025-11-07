library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

####SET PARAMETERS####

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#define location to save QC output data
input_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/1_platemaps_resources"
save_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/qc_plots"

#guide position for neg controls on locus plot
negguidepos = 22150000
locguidepos = 22160000
celltype = "smc"


#define the replicate numbers, metadata, etc to be assigned to each plate
dxl04_repnum1 = 15
dxl04_repnum2 = 16
dxl04_seqnum = 8

dxl05_repnum1 = 17
dxl05_repnum2 = 18
dxl05_seqnum = 9

dxl02_repnum1 = 19
dxl02_repnum2 = 20
dxl02_seqnum = 10

dxl03_repnum1 = 21
dxl03_repnum2 = 22
dxl03_seqnum = 11

dxl01_repnum1 = 23
dxl01_repnum2 = 24
dxl01_seqnum = 12

dxk99_repnum1 = 25
dxk99_repnum2 = 26
dxk99_seqnum = 13

dxk98_repnum1 = 27
dxk98_repnum2 = 28
dxk98_seqnum = 14

dxk97_repnum1 = 29
dxk97_repnum2 = 30
dxk97_seqnum = 15

#read in barcode keys to associate read data and metadata

dxl04_barcode_key <- read.csv(paste0(input_path,"/dxl04_key.csv"))
dxl05_barcode_key <- read.csv(paste0(input_path,"/dxl05_key.csv"))
dxl02_barcode_key <- read.csv(paste0(input_path,"/dxl02_key.csv"))
dxl03_barcode_key <- read.csv(paste0(input_path,"/dxl03_key.csv"))
dxl01_barcode_key <- read.csv(paste0(input_path,"/dxl01_key.csv"))
dxk99_barcode_key <- read.csv(paste0(input_path,"/dxk99_key.csv"))
dxk98_barcode_key <- read.csv(paste0(input_path,"/dxk98_key.csv"))
dxk97_barcode_key <- read.csv(paste0(input_path,"/dxk97_key.csv"))

#read in the screen key, which includes all screen and control guides
#screen_key <- read.csv("platemaps/9p21_Screen_sgRNA_key.csv")
screen_key <- read.csv(paste0(input_path,"/9p21_Screen_sgRNA_key.csv"))

#merge each barcode key with additional metadata values from the full screen
#these include guide name and position, etc
dxl04_barcode_key <- merge(dxl04_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxl04_barcode_key <- dxl04_barcode_key[!duplicated(dxl04_barcode_key$well),]

dxl05_barcode_key <- merge(dxl05_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxl05_barcode_key <- dxl05_barcode_key[!duplicated(dxl05_barcode_key$well),]

dxl02_barcode_key <- merge(dxl02_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxl02_barcode_key <- dxl02_barcode_key[!duplicated(dxl02_barcode_key$well),]

dxl03_barcode_key <- merge(dxl03_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxl03_barcode_key <- dxl03_barcode_key[!duplicated(dxl03_barcode_key$well),]

dxl01_barcode_key <- merge(dxl01_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxl01_barcode_key <- dxl01_barcode_key[!duplicated(dxl01_barcode_key$well),]

dxk99_barcode_key <- merge(dxk99_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxk99_barcode_key <- dxk99_barcode_key[!duplicated(dxk99_barcode_key$well),]

dxk98_barcode_key <- merge(dxk98_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxk98_barcode_key <- dxk98_barcode_key[!duplicated(dxk98_barcode_key$well),]

dxk97_barcode_key <- merge(dxk97_barcode_key, screen_key[,c("guide_id", "GuidePosition", "GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
dxk97_barcode_key <- dxk97_barcode_key[!duplicated(dxk97_barcode_key$well),]

#for each plate, add a metadata factor for guide ID; screen, locus control, or scrambled control
#this will allow effective filtering of the screen guides later on
dxl04_barcode_key <- dxl04_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxl05_barcode_key <- dxl05_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxl02_barcode_key <- dxl02_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxl03_barcode_key <- dxl03_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxl01_barcode_key <- dxl01_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxk99_barcode_key <- dxk99_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxk98_barcode_key <- dxk98_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    #ifelse(is.na(GuidePosition), "Selection control",
                                    ifelse(sequence == "", "Selection control",
                                           "Screen")))
  )

dxk97_barcode_key <- dxk97_barcode_key %>%
  mutate(cellType = celltype,
         guide_ctrl = ifelse(guide == "Neg-9p21 regions", "Locus control",
                             ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                    #ifelse(is.na(GuidePosition), "Selection control",
                                    ifelse(sequence == "", "Selection control",
                                                  "Screen")))
  )

#import demultiplexed STARsolo data using Seurat's ReadMtx function
#this is redundant with original QC files (but better formatted here)
dxl04_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_1/dxl04_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_1/dxl04_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_1/dxl04_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl04_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_2/dxl04_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_2/dxl04_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/STARsolo_2/dxl04_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl05_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_1/dxl05_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_1/dxl05_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_1/dxl05_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl05_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_2/dxl05_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_2/dxl05_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/STARsolo_2/dxl05_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl03_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_1/dxl03_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_1/dxl03_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_1/dxl03_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl03_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_2/dxl03_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_2/dxl03_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/STARsolo_2/dxl03_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl02_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_1/dxl02_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_1/dxl02_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_1/dxl02_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl02_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_2/dxl02_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_2/dxl02_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/STARsolo_2/dxl02_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl01_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_1/dxl01_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_1/dxl01_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_1/dxl01_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxl01_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_2/dxl01_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_2/dxl01_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/STARsolo_2/dxl01_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk99_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_1/dxk99_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_1/dxk99_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_1/dxk99_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk99_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_2/dxk99_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_2/dxk99_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/STARsolo_2/dxk99_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk98_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_1/dxk98_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_1/dxk98_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_1/dxk98_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk98_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_2/dxk98_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_2/dxk98_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/STARsolo_2/dxk98_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk97_1_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_1/dxk97_smc_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_1/dxk97_smc_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_1/dxk97_smc_1Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)

dxk97_2_smc.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_2/dxk97_smc_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_2/dxk97_smc_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/STARsolo_2/dxk97_smc_2Solo.out/Gene/filtered/features.tsv",
                            cell.column = 1,
                            feature.column = 2,
                            cell.sep = "\t",
                            feature.sep = "\t",
                            skip.cell = 0,
                            skip.feature = 0,
                            mtx.transpose = FALSE,
                            unique.features = TRUE,
                            strip.suffix = FALSE
)


#convert genecount data to a dataframe by way of R matrix
dxl04_genecounts1 <- as.matrix(dxl04_1_smc.data)
dxl04_genecounts2 <- as.matrix(dxl04_2_smc.data)
dxl04_genecounts1 <- as.data.frame(dxl04_genecounts1)
dxl04_genecounts2 <- as.data.frame(dxl04_genecounts2)

dxl05_genecounts1 <- as.matrix(dxl05_1_smc.data)
dxl05_genecounts2 <- as.matrix(dxl05_2_smc.data)
dxl05_genecounts1 <- as.data.frame(dxl05_genecounts1)
dxl05_genecounts2 <- as.data.frame(dxl05_genecounts2)

dxl03_genecounts1 <- as.matrix(dxl03_1_smc.data)
dxl03_genecounts2 <- as.matrix(dxl03_2_smc.data)
dxl03_genecounts1 <- as.data.frame(dxl03_genecounts1)
dxl03_genecounts2 <- as.data.frame(dxl03_genecounts2)

dxl02_genecounts1 <- as.matrix(dxl02_1_smc.data)
dxl02_genecounts2 <- as.matrix(dxl02_2_smc.data)
dxl02_genecounts1 <- as.data.frame(dxl02_genecounts1)
dxl02_genecounts2 <- as.data.frame(dxl02_genecounts2)

dxl01_genecounts1 <- as.matrix(dxl01_1_smc.data)
dxl01_genecounts2 <- as.matrix(dxl01_2_smc.data)
dxl01_genecounts1 <- as.data.frame(dxl01_genecounts1)
dxl01_genecounts2 <- as.data.frame(dxl01_genecounts2)

dxk99_genecounts1 <- as.matrix(dxk99_1_smc.data)
dxk99_genecounts2 <- as.matrix(dxk99_2_smc.data)
dxk99_genecounts1 <- as.data.frame(dxk99_genecounts1)
dxk99_genecounts2 <- as.data.frame(dxk99_genecounts2)

dxk98_genecounts1 <- as.matrix(dxk98_1_smc.data)
dxk98_genecounts2 <- as.matrix(dxk98_2_smc.data)
dxk98_genecounts1 <- as.data.frame(dxk98_genecounts1)
dxk98_genecounts2 <- as.data.frame(dxk98_genecounts2)

dxk97_genecounts1 <- as.matrix(dxk97_1_smc.data)
dxk97_genecounts2 <- as.matrix(dxk97_2_smc.data)
dxk97_genecounts1 <- as.data.frame(dxk97_genecounts1)
dxk97_genecounts2 <- as.data.frame(dxk97_genecounts2)


####QC: extract depth, %ribo, %mito####

#extract read counts, %mito, and %ribo per well, and export as csv
#repeat for each plate

#get readcounts for dxl04
readcounts <- data.frame(depth = colSums(dxl04_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxl04_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxl04_repnum1,
         seq_id = dxl04_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxl04_repnum2,
         seq_id = dxl04_seqnum,
  )

mt1 <- data.frame(dxl04_genecounts1[grep('^MT-',rownames(dxl04_genecounts1)),])
mt2 <- data.frame(dxl04_genecounts2[grep('^MT-',rownames(dxl04_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxl04_genecounts1[grep('^RPS|^RPL',rownames(dxl04_genecounts1)),])
ribo2 <- data.frame(dxl04_genecounts2[grep('^RPS|^RPL',rownames(dxl04_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl04_readcounts <- merge(dxl04_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl04_readcounts2 <- merge(dxl04_barcode_key, readcounts2, by="bc", all=T)

#write.csv(readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/dxl04_smc_qc_gene_counts.csv")

#get readcounts and metadata for dxl05
readcounts <- data.frame(depth = colSums(dxl05_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxl05_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxl05_repnum1,
         seq_id = dxl05_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxl05_repnum2,
         seq_id = dxl05_seqnum,
  )

mt1 <- data.frame(dxl05_genecounts1[grep('^MT-',rownames(dxl05_genecounts1)),])
mt2 <- data.frame(dxl05_genecounts2[grep('^MT-',rownames(dxl05_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxl05_genecounts1[grep('^RPS|^RPL',rownames(dxl05_genecounts1)),])
ribo2 <- data.frame(dxl05_genecounts2[grep('^RPS|^RPL',rownames(dxl05_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl05_readcounts <- merge(dxl05_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl05_readcounts2 <- merge(dxl05_barcode_key, readcounts2, by="bc", all=T)

#write.csv(readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/dxl05_smc_qc_gene_counts.csv")

#get readcounts and metadata for dxl03
readcounts <- data.frame(depth = colSums(dxl03_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxl03_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxl03_repnum1,
         seq_id = dxl03_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxl03_repnum2,
         seq_id = dxl03_seqnum,
  )

mt1 <- data.frame(dxl03_genecounts1[grep('^MT-',rownames(dxl03_genecounts1)),])
mt2 <- data.frame(dxl03_genecounts2[grep('^MT-',rownames(dxl03_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxl03_genecounts1[grep('^RPS|^RPL',rownames(dxl03_genecounts1)),])
ribo2 <- data.frame(dxl03_genecounts2[grep('^RPS|^RPL',rownames(dxl03_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl03_readcounts <- merge(dxl03_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl03_readcounts2 <- merge(dxl03_barcode_key, readcounts2, by="bc", all=T)

#write.csv(readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/dxl03_smc_qc_gene_counts.csv")

#get readcounts and metadata for dxl02
readcounts <- data.frame(depth = colSums(dxl02_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxl02_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxl02_repnum1,
         seq_id = dxl02_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxl02_repnum2,
         seq_id = dxl02_seqnum,
  )

mt1 <- data.frame(dxl02_genecounts1[grep('^MT-',rownames(dxl02_genecounts1)),])
mt2 <- data.frame(dxl02_genecounts2[grep('^MT-',rownames(dxl02_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxl02_genecounts1[grep('^RPS|^RPL',rownames(dxl02_genecounts1)),])
ribo2 <- data.frame(dxl02_genecounts2[grep('^RPS|^RPL',rownames(dxl02_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl02_readcounts <- merge(dxl02_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl02_readcounts2 <- merge(dxl02_barcode_key, readcounts2, by="bc", all=T)

#get readcounts and metadata  for dxl01
readcounts <- data.frame(depth = colSums(dxl01_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxl01_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxl01_repnum1,
         seq_id = dxl01_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxl01_repnum2,
         seq_id = dxl01_seqnum,
  )

mt1 <- data.frame(dxl01_genecounts1[grep('^MT-',rownames(dxl01_genecounts1)),])
mt2 <- data.frame(dxl01_genecounts2[grep('^MT-',rownames(dxl01_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxl01_genecounts1[grep('^RPS|^RPL',rownames(dxl01_genecounts1)),])
ribo2 <- data.frame(dxl01_genecounts2[grep('^RPS|^RPL',rownames(dxl01_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl01_readcounts <- merge(dxl01_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxl01_readcounts2 <- merge(dxl01_barcode_key, readcounts2, by="bc", all=T)

#get readcounts and metadata for dxk99
readcounts <- data.frame(depth = colSums(dxk99_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxk99_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxk99_repnum1,
         seq_id = dxk99_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxk99_repnum2,
         seq_id = dxk99_seqnum,
  )

mt1 <- data.frame(dxk99_genecounts1[grep('^MT-',rownames(dxk99_genecounts1)),])
mt2 <- data.frame(dxk99_genecounts2[grep('^MT-',rownames(dxk99_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxk99_genecounts1[grep('^RPS|^RPL',rownames(dxk99_genecounts1)),])
ribo2 <- data.frame(dxk99_genecounts2[grep('^RPS|^RPL',rownames(dxk99_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk99_readcounts <- merge(dxk99_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk99_readcounts2 <- merge(dxk99_barcode_key, readcounts2, by="bc", all=T)

#get readcounts and metadata dxk98
readcounts <- data.frame(depth = colSums(dxk98_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxk98_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxk98_repnum1,
         seq_id = dxk98_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxk98_repnum2,
         seq_id = dxk98_seqnum,
  )

mt1 <- data.frame(dxk98_genecounts1[grep('^MT-',rownames(dxk98_genecounts1)),])
mt2 <- data.frame(dxk98_genecounts2[grep('^MT-',rownames(dxk98_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxk98_genecounts1[grep('^RPS|^RPL',rownames(dxk98_genecounts1)),])
ribo2 <- data.frame(dxk98_genecounts2[grep('^RPS|^RPL',rownames(dxk98_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk98_readcounts <- merge(dxk98_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk98_readcounts2 <- merge(dxk98_barcode_key, readcounts2, by="bc", all=T)

#get readcounts and metadata for dxk97
readcounts <- data.frame(depth = colSums(dxk97_genecounts1))
readcounts2 <- data.frame(depth = colSums(dxk97_genecounts2))

readcounts <- readcounts %>%
  mutate(rep_id = dxk97_repnum1,
         seq_id = dxk97_seqnum,
  )

readcounts2 <- readcounts2 %>%
  mutate(rep_id = dxk97_repnum2,
         seq_id = dxk97_seqnum,
  )

mt1 <- data.frame(dxk97_genecounts1[grep('^MT-',rownames(dxk97_genecounts1)),])
mt2 <- data.frame(dxk97_genecounts2[grep('^MT-',rownames(dxk97_genecounts2)),])
mt_counts <- data.frame(mt_count = colSums(mt1))
mt_counts2 <- data.frame(mt_count = colSums(mt2))

ribo1 <- data.frame(dxk97_genecounts1[grep('^RPS|^RPL',rownames(dxk97_genecounts1)),])
ribo2 <- data.frame(dxk97_genecounts2[grep('^RPS|^RPL',rownames(dxk97_genecounts2)),])
ribo_counts <- data.frame(ribo_count = colSums(ribo1))
ribo_counts2 <- data.frame(ribo_count = colSums(ribo1))

readcounts <- merge(readcounts, mt_counts, by="row.names")
colnames(readcounts)[1] <- "bc"
readcounts <- merge(readcounts, ribo_counts, by.x = "bc", by.y="row.names")
readcounts <- readcounts %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk97_readcounts <- merge(dxk97_barcode_key, readcounts, by="bc", all=T)

readcounts2 <- merge(readcounts2, mt_counts2, by="row.names")
colnames(readcounts2)[1] <- "bc"
readcounts2 <- merge(readcounts2, ribo_counts2, by.x = "bc", by.y="row.names")
readcounts2 <- readcounts2 %>% mutate(
  mtper = depth/mt_count,
  riboper = depth/ribo_count,
)
dxk97_readcounts2 <- merge(dxk97_barcode_key, readcounts2, by="bc", all=T)

#write.csv(readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/dxl02_smc_qc_gene_counts.csv")

####collect counts of 9p21 genes####

#incorporate metadata into counts for individual and pooled runs
#this function pulls out counts of one gene and appends it to the rea counts file
#this will be run for each 9p21 gene
getcounts_rep <- function (a,b,d) {
  counts1 <- t(subset(a, rownames(a) %in% c(b)))
  d <- merge(d, counts1, by.x="bc", by.y="row.names", all=T)
}

dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "CDKN2B-AS1", dxl04_readcounts)
dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "MTAP", dxl04_readcounts)
dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "CDKN2A", dxl04_readcounts)
dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "CDKN2B", dxl04_readcounts)
dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "DMRTA1", dxl04_readcounts)
dxl04_readcounts <- getcounts_rep(dxl04_genecounts1, "MIR31HG", dxl04_readcounts)

dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "CDKN2B-AS1", dxl04_readcounts2)
dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "MTAP", dxl04_readcounts2)
dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "CDKN2A", dxl04_readcounts2)
dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "CDKN2B", dxl04_readcounts2)
dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "DMRTA1", dxl04_readcounts2)
dxl04_readcounts2 <- getcounts_rep(dxl04_genecounts2, "MIR31HG", dxl04_readcounts2)

write.csv(dxl04_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/dxl04_smc_rep1_raw.csv")
write.csv(dxl04_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/dxl04_smc_rep2_raw.csv")

dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "CDKN2B-AS1", dxl05_readcounts)
dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "MTAP", dxl05_readcounts)
dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "CDKN2A", dxl05_readcounts)
dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "CDKN2B", dxl05_readcounts)
dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "DMRTA1", dxl05_readcounts)
dxl05_readcounts <- getcounts_rep(dxl05_genecounts1, "MIR31HG", dxl05_readcounts)

dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "CDKN2B-AS1", dxl05_readcounts2)
dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "MTAP", dxl05_readcounts2)
dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "CDKN2A", dxl05_readcounts2)
dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "CDKN2B", dxl05_readcounts2)
dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "DMRTA1", dxl05_readcounts2)
dxl05_readcounts2 <- getcounts_rep(dxl05_genecounts2, "MIR31HG", dxl05_readcounts2)

write.csv(dxl05_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/dxl05_smc_rep1_raw.csv")
write.csv(dxl05_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/dxl05_smc_rep2_raw.csv")

dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "CDKN2B-AS1", dxl03_readcounts)
dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "MTAP", dxl03_readcounts)
dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "CDKN2A", dxl03_readcounts)
dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "CDKN2B", dxl03_readcounts)
dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "DMRTA1", dxl03_readcounts)
dxl03_readcounts <- getcounts_rep(dxl03_genecounts1, "MIR31HG", dxl03_readcounts)

dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "CDKN2B-AS1", dxl03_readcounts2)
dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "MTAP", dxl03_readcounts2)
dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "CDKN2A", dxl03_readcounts2)
dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "CDKN2B", dxl03_readcounts2)
dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "DMRTA1", dxl03_readcounts2)
dxl03_readcounts2 <- getcounts_rep(dxl03_genecounts2, "MIR31HG", dxl03_readcounts2)

write.csv(dxl03_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/dxl03_smc_rep1_raw.csv")
write.csv(dxl03_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/dxl03_smc_rep2_raw.csv")

dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "CDKN2B-AS1", dxl02_readcounts)
dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "MTAP", dxl02_readcounts)
dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "CDKN2A", dxl02_readcounts)
dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "CDKN2B", dxl02_readcounts)
dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "DMRTA1", dxl02_readcounts)
dxl02_readcounts <- getcounts_rep(dxl02_genecounts1, "MIR31HG", dxl02_readcounts)

dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "CDKN2B-AS1", dxl02_readcounts2)
dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "MTAP", dxl02_readcounts2)
dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "CDKN2A", dxl02_readcounts2)
dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "CDKN2B", dxl02_readcounts2)
dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "DMRTA1", dxl02_readcounts2)
dxl02_readcounts2 <- getcounts_rep(dxl02_genecounts2, "MIR31HG", dxl02_readcounts2)

write.csv(dxl02_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/dxl02_smc_rep1_raw.csv")
write.csv(dxl02_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/dxl02_smc_rep2_raw.csv")

dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "CDKN2B-AS1", dxl01_readcounts)
dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "MTAP", dxl01_readcounts)
dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "CDKN2A", dxl01_readcounts)
dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "CDKN2B", dxl01_readcounts)
dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "DMRTA1", dxl01_readcounts)
dxl01_readcounts <- getcounts_rep(dxl01_genecounts1, "MIR31HG", dxl01_readcounts)

dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "CDKN2B-AS1", dxl01_readcounts2)
dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "MTAP", dxl01_readcounts2)
dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "CDKN2A", dxl01_readcounts2)
dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "CDKN2B", dxl01_readcounts2)
dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "DMRTA1", dxl01_readcounts2)
dxl01_readcounts2 <- getcounts_rep(dxl01_genecounts2, "MIR31HG", dxl01_readcounts2)

write.csv(dxl01_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/dxl01_smc_rep1_raw.csv")
write.csv(dxl01_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/dxl01_smc_rep2_raw.csv")

dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "CDKN2B-AS1", dxk99_readcounts)
dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "MTAP", dxk99_readcounts)
dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "CDKN2A", dxk99_readcounts)
dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "CDKN2B", dxk99_readcounts)
dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "DMRTA1", dxk99_readcounts)
dxk99_readcounts <- getcounts_rep(dxk99_genecounts1, "MIR31HG", dxk99_readcounts)

dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "CDKN2B-AS1", dxk99_readcounts2)
dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "MTAP", dxk99_readcounts2)
dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "CDKN2A", dxk99_readcounts2)
dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "CDKN2B", dxk99_readcounts2)
dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "DMRTA1", dxk99_readcounts2)
dxk99_readcounts2 <- getcounts_rep(dxk99_genecounts1, "MIR31HG", dxk99_readcounts2)

write.csv(dxk99_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/dxk99_smc_rep1_raw.csv")
write.csv(dxk99_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/dxk99_smc_rep2_raw.csv")

dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "CDKN2B-AS1", dxk98_readcounts)
dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "MTAP", dxk98_readcounts)
dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "CDKN2A", dxk98_readcounts)
dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "CDKN2B", dxk98_readcounts)
dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "DMRTA1", dxk98_readcounts)
dxk98_readcounts <- getcounts_rep(dxk98_genecounts1, "MIR31HG", dxk98_readcounts)

dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "CDKN2B-AS1", dxk98_readcounts2)
dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "MTAP", dxk98_readcounts2)
dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "CDKN2A", dxk98_readcounts2)
dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "CDKN2B", dxk98_readcounts2)
dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "DMRTA1", dxk98_readcounts2)
dxk98_readcounts2 <- getcounts_rep(dxk98_genecounts1, "MIR31HG", dxk98_readcounts2)

write.csv(dxk98_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/dxk98_smc_rep1_raw.csv")
write.csv(dxk98_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/dxk98_smc_rep2_raw.csv")

dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "CDKN2B-AS1", dxk97_readcounts)
dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "MTAP", dxk97_readcounts)
dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "CDKN2A", dxk97_readcounts)
dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "CDKN2B", dxk97_readcounts)
dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "DMRTA1", dxk97_readcounts)
dxk97_readcounts <- getcounts_rep(dxk97_genecounts1, "MIR31HG", dxk97_readcounts)

dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "CDKN2B-AS1", dxk97_readcounts2)
dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "MTAP", dxk97_readcounts2)
dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "CDKN2A", dxk97_readcounts2)
dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "CDKN2B", dxk97_readcounts2)
dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "DMRTA1", dxk97_readcounts2)
dxk97_readcounts2 <- getcounts_rep(dxk97_genecounts1, "MIR31HG", dxk97_readcounts2)

write.csv(dxk97_readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/dxk97_smc_rep1_raw.csv")
write.csv(dxk97_readcounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/dxk97_smc_rep2_raw.csv")

####TPM normalization: export TPM counts in separate document####

#a=destination file to convert
#b=file with metadata and counts
tpm_convert <- function(a,b) {
  a <- b
  a %>% mutate(
    `CDKN2B-AS1` = (10^6*`CDKN2B-AS1`)/depth,
    MTAP = (10^6*MTAP)/depth,
    CDKN2A = (10^6*CDKN2A)/depth,
    CDKN2B = (10^6*CDKN2B)/depth,
    DMRTA1 = (10^6*DMRTA1)/depth,
    MIR31HG = (10^6*MIR31HG)/depth
  )
  
}

dxl04_readcounts_tpm <- tpm_convert(dxl04_readcounts_tpm,dxl04_readcounts)
dxl04_readcounts2_tpm <- tpm_convert(dxl04_readcounts2_tpm,dxl04_readcounts2)

write.csv(dxl04_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/dxl04_smc_rep1_tpm.csv")
write.csv(dxl04_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/dxl04_smc_rep2_tpm.csv")

dxl05_readcounts_tpm <- tpm_convert(dxl05_readcounts_tpm,dxl05_readcounts)
dxl05_readcounts2_tpm <- tpm_convert(dxl05_readcounts2_tpm,dxl05_readcounts2)

write.csv(dxl05_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/dxl05_smc_rep1_tpm.csv")
write.csv(dxl05_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/dxl05_smc_rep2_tpm.csv")

dxl03_readcounts_tpm <- tpm_convert(dxl03_readcounts_tpm,dxl03_readcounts)
dxl03_readcounts2_tpm <- tpm_convert(dxl03_readcounts2_tpm,dxl03_readcounts2)

write.csv(dxl03_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/dxl03_smc_rep1_tpm.csv")
write.csv(dxl03_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/dxl03_smc_rep2_tpm.csv")

dxl02_readcounts_tpm <- tpm_convert(dxl02_readcounts_tpm,dxl02_readcounts)
dxl02_readcounts2_tpm <- tpm_convert(dxl02_readcounts2_tpm,dxl02_readcounts2)

write.csv(dxl02_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/dxl02_smc_rep1_tpm.csv")
write.csv(dxl02_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/dxl02_smc_rep2_tpm.csv")

dxl01_readcounts_tpm <- tpm_convert(dxl01_readcounts_tpm,dxl01_readcounts)
dxl01_readcounts2_tpm <- tpm_convert(dxl01_readcounts2_tpm,dxl01_readcounts2)

write.csv(dxl01_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/dxl01_smc_rep1_tpm.csv")
write.csv(dxl01_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/dxl01_smc_rep2_tpm.csv")

dxk99_readcounts_tpm <- tpm_convert(dxk99_readcounts_tpm,dxk99_readcounts)
dxk99_readcounts2_tpm <- tpm_convert(dxk99_readcounts2_tpm,dxk99_readcounts2)

write.csv(dxk99_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/dxk99_smc_rep1_tpm.csv")
write.csv(dxk99_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/dxk99_smc_rep2_tpm.csv")

dxk98_readcounts_tpm <- tpm_convert(dxk98_readcounts_tpm,dxk98_readcounts)
dxk98_readcounts2_tpm <- tpm_convert(dxk98_readcounts2_tpm,dxk98_readcounts2)

write.csv(dxk98_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/dxk98_smc_rep1_tpm.csv")
write.csv(dxk98_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/dxk98_smc_rep2_tpm.csv")

dxk97_readcounts_tpm <- tpm_convert(dxk97_readcounts_tpm,dxk97_readcounts)
dxk97_readcounts2_tpm <- tpm_convert(dxk97_readcounts2_tpm,dxk97_readcounts2)

write.csv(dxk97_readcounts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/dxk97_smc_rep1_tpm.csv")
write.csv(dxk97_readcounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/dxk97_smc_rep2_tpm.csv")


####convert bc label to guide id label in whole-transcriptome matrix####
#this means all gene count matrices can be combined for RNA-Seq
#for dxl04

#sort the barcode key alphabetically by barcode sequence (as it appears in the gene count matrix)
dxl04_barcode_key <- dxl04_barcode_key %>% arrange(bc)
#extract all present barcodes in alphabetical order from the gene count matrix
test_cols1 <- data.frame(colnames(dxl04_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxl04_genecounts2))
colnames(test_cols2) <- c("bc")

#create a data frame that combines the barcodes (in order) and the ID column for individual guides
test_cols1 <- merge(test_cols1,dxl04_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl04_genecounts1) <- test_cols1$guide_id
#repeat the above for the second replicate
test_cols2 <- merge(test_cols2,dxl04_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl04_genecounts2) <- test_cols2$guide_id
colnames(dxl04_genecounts1) <- paste0(colnames(dxl04_genecounts1),"_rep1")
colnames(dxl04_genecounts2) <- paste0(colnames(dxl04_genecounts2),"_rep2")

write.csv(dxl04_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/raw_counts_dxl04_smc1.csv")
write.csv(dxl04_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl04/raw_counts_dxl04_smc2.csv")

#for dxl05
dxl05_barcode_key <- dxl05_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxl05_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxl05_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxl05_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl05_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxl05_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl05_genecounts2) <- test_cols2$guide_id
colnames(dxl05_genecounts1) <- paste0(colnames(dxl05_genecounts1),"_rep1")
colnames(dxl05_genecounts2) <- paste0(colnames(dxl05_genecounts2),"_rep2")

write.csv(dxl05_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/raw_counts_dxl05_smc1.csv")
write.csv(dxl05_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl05/raw_counts_dxl05_smc2.csv")

#for dxl03
dxl03_barcode_key <- dxl03_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxl03_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxl03_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxl03_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl03_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxl03_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl03_genecounts2) <- test_cols2$guide_id
colnames(dxl03_genecounts1) <- paste0(colnames(dxl03_genecounts1),"_rep1")
colnames(dxl03_genecounts2) <- paste0(colnames(dxl03_genecounts2),"_rep2")

write.csv(dxl03_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/raw_counts_dxl03_smc1.csv")
write.csv(dxl03_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl03/raw_counts_dxl03_smc2.csv")

#for dxl02
dxl02_barcode_key <- dxl02_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxl02_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxl02_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxl02_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl02_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxl02_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl02_genecounts2) <- test_cols2$guide_id
colnames(dxl02_genecounts1) <- paste0(colnames(dxl02_genecounts1),"_rep1")
colnames(dxl02_genecounts2) <- paste0(colnames(dxl02_genecounts2),"_rep2")

write.csv(dxl02_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/raw_counts_dxl02_smc1.csv")
write.csv(dxl02_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl02/raw_counts_dxl02_smc2.csv")

#for dxl01
dxl01_barcode_key <- dxl01_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxl01_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxl01_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxl01_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl01_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxl01_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxl01_genecounts2) <- test_cols2$guide_id
colnames(dxl01_genecounts1) <- paste0(colnames(dxl01_genecounts1),"_rep1")
colnames(dxl01_genecounts2) <- paste0(colnames(dxl01_genecounts2),"_rep2")

write.csv(dxl01_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/raw_counts_dxl01_smc1.csv")
write.csv(dxl01_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxl01/raw_counts_dxl01_smc2.csv")

#for dxk99
dxk99_barcode_key <- dxk99_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxk99_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxk99_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxk99_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk99_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxk99_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk99_genecounts2) <- test_cols2$guide_id
colnames(dxk99_genecounts1) <- paste0(colnames(dxk99_genecounts2),"_rep1")
colnames(dxk99_genecounts2) <- paste0(colnames(dxk99_genecounts2),"_rep2")

write.csv(dxk99_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/raw_counts_dxk99_smc1.csv")
write.csv(dxk99_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk99/raw_counts_dxk99_smc2.csv")

#for dxk98
dxk98_barcode_key <- dxk98_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxk98_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxk98_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxk98_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk98_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxk98_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk98_genecounts2) <- test_cols2$guide_id
colnames(dxk98_genecounts1) <- paste0(colnames(dxk98_genecounts1),"_rep1")
colnames(dxk98_genecounts2) <- paste0(colnames(dxk98_genecounts2),"_rep2")

write.csv(dxk98_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/raw_counts_dxk98_smc1.csv")
write.csv(dxk98_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk98/raw_counts_dxk98_smc2.csv")

#for dxk97
dxk97_barcode_key <- dxk97_barcode_key %>% arrange(bc)
test_cols1 <- data.frame(colnames(dxk97_genecounts1))
colnames(test_cols1) <- c("bc")
test_cols2 <- data.frame(colnames(dxk97_genecounts2))
colnames(test_cols2) <- c("bc")

test_cols1 <- merge(test_cols1,dxk97_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk97_genecounts1) <- test_cols1$guide_id
test_cols2 <- merge(test_cols2,dxk97_barcode_key[,c("bc", "guide_id")], by.x="bc", by.y="bc", all.x =T)
colnames(dxk97_genecounts2) <- test_cols2$guide_id
colnames(dxk97_genecounts1) <- paste0(colnames(dxk97_genecounts1),"_rep1")
colnames(dxk97_genecounts2) <- paste0(colnames(dxk97_genecounts2),"_rep2")

write.csv(dxk97_genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/raw_counts_dxk97_smc1.csv")
write.csv(dxk97_genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/smc_dxk97/raw_counts_dxk97_smc2.csv")






