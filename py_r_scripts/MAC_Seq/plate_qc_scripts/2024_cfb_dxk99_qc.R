library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


####SET PARAMETERS####

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#save_path = "output/cfb_dxk99"
save_path = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/qc_plots"

#guide position for neg controls on locus plot
negguidepos = 22150000
locguidepos = 22160000
celltype = "cfb"
repnum1 = 21
repnum2 = 22
seqnum = 11

barcode_key <- read.csv("raw_data_now_on_helium/raw_data/dxk99_key.csv")
screen_key <- read.csv("platemaps/9p21_Screen_sgRNA_key.csv")

barcode_key <- merge(barcode_key, screen_key[,c("GuidePosition","GuideSeq")], by.x="sequence", by.y="GuideSeq", all.x =T)
barcode_key <- barcode_key[!duplicated(barcode_key$well),]

barcode_key <- barcode_key %>%
  mutate(cellType = celltype,
         GuidePosition = ifelse(guide == "NEG_CONTROL", negguidepos,
                                ifelse(guide == "Neg-9p21 regions", locguidepos,
                                       GuidePosition)),
         guide_ctrl = ifelse(is.na(GuidePosition), "Selection control",
                             ifelse(guide == "Neg-9p21 regions", "Locus control",
                                    ifelse(guide == "NEG_CONTROL", "Scrambled control",
                                           "Screen")))
  )

####READ_MATRICES####

#import STARsolo data 
dxk99_1_cfb.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_1/dxk99_cfb_1Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_1/dxk99_cfb_1Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_1/dxk99_cfb_1Solo.out/Gene/filtered/features.tsv",
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

dxk99_2_cfb.data <- ReadMtx(mtx = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_2/dxk99_cfb_2Solo.out/Gene/filtered/matrix.mtx",
                            cells = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_2/dxk99_cfb_2Solo.out/Gene/filtered/barcodes.tsv", 
                            features = "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/STARsolo_2/dxk99_cfb_2Solo.out/Gene/filtered/features.tsv",
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

genecounts1 <- as.matrix(dxk99_1_cfb.data)
genecounts2 <- as.matrix(dxk99_2_cfb.data)

write.csv(genecounts1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/raw_counts_dxk99_cfb1.csv")

write.csv(genecounts2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/raw_counts_dxk99_cfb2.csv")

genecounts1 <- as.data.frame(genecounts1)
genecounts2 <- as.data.frame(genecounts2)

###collect counts of 9p21 genes###

#define equation for linear regression plotting
lm_eq <- function(z) {substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                 
                                 list(        a = as.character(as.data.frame(format(coef(z)[1], digits = 3))),
                                              
                                              b = as.character(as.data.frame(format(coef(z)[2], digits = 3))),
                                              
                                              r2 = format(summary(z)$r.squared, digits = 3)))
}

#incorporate metadata into counts for individual and pooled runs

getcounts <- function (a,b,d,e) {
  counts1 <- t(subset(a, rownames(a) %in% c(d)))
  counts2 <- t(subset(b, rownames(b) %in% c(d)))
  e <- merge(e, counts1, by.x="bc", by.y="row.names", all=T)
  e <- merge(e, counts2, by.x="bc", by.y="row.names", all=T)
}

getcounts_rep <- function (a,b,d) {
  counts1 <- t(subset(a, rownames(a) %in% c(b)))
  d <- merge(d, counts1, by.x="bc", by.y="row.names", all=T)
}

counts <- barcode_key

counts_rep1 <- barcode_key %>%
  mutate(rep_id = repnum1,
         seq_id = seqnum,
  )
counts_rep1 <- merge(counts_rep1, colSums(genecounts1), by.x="bc",by.y="row.names", all=T)

counts_rep2 <- barcode_key %>%
  mutate(rep_id = repnum2,
         seq_id = seqnum,
  )
counts_rep2 <- merge(counts_rep2, colSums(genecounts2), by.x="bc",by.y="row.names", all=T)

#get counts of all 9p21 genes and append to sample info
counts <- getcounts(genecounts1,genecounts2,"CDKN2B-AS1",counts)
counts <- getcounts(genecounts1,genecounts2,"MTAP",counts)
counts <- getcounts(genecounts1,genecounts2,"CDKN2A",counts)
counts <- getcounts(genecounts1,genecounts2,"CDKN2B",counts)
counts <- getcounts(genecounts1,genecounts2,"DMRTA1",counts)
counts <- getcounts(genecounts1,genecounts2,"MIR31HG",counts)
#repeat for first plate only
counts_rep1 <- getcounts_rep(genecounts1, "CDKN2B-AS1", counts_rep1)
counts_rep1 <- getcounts_rep(genecounts1, "MTAP", counts_rep1)
counts_rep1 <- getcounts_rep(genecounts1, "CDKN2A", counts_rep1)
counts_rep1 <- getcounts_rep(genecounts1, "CDKN2B", counts_rep1)
counts_rep1 <- getcounts_rep(genecounts1, "DMRTA1", counts_rep1)
counts_rep1 <- getcounts_rep(genecounts1, "MIR31HG", counts_rep1)
#and for second plate only
counts_rep2 <- getcounts_rep(genecounts2, "CDKN2B-AS1", counts_rep2)
counts_rep2 <- getcounts_rep(genecounts2, "MTAP", counts_rep2)
counts_rep2 <- getcounts_rep(genecounts2, "CDKN2A", counts_rep2)
counts_rep2 <- getcounts_rep(genecounts2, "CDKN2B", counts_rep2)
counts_rep2 <- getcounts_rep(genecounts2, "DMRTA1", counts_rep2)
counts_rep2 <- getcounts_rep(genecounts2, "MIR31HG", counts_rep2)

write.csv(counts_rep1, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_rep1_raw.csv")
write.csv(counts_rep2, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_rep2_raw.csv")

####TPM_NORMALIZATION####

#genesums <- as.data.frame(colSums(dxl04_1.data))
#colnames(genesums) <- c("counts")

tpm <- function(x) {10^6*x/colSums(as.matrix(x))}
logtpm <- function(y) {log10(1+(10^6*y/colSums(as.matrix(y))))}

#calculate TPM values for all 9p21 genes
genecounts1_tpm <- genecounts1 %>% mutate_each(tpm)
genecounts2_tpm <- genecounts2 %>% mutate_each(tpm)

write.csv(genecounts1_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/tpm_dxk99_cfb1_all.csv")
write.csv(genecounts2_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/tpm_dxk99_cfb2_all.csv")

#define a function to label 9p21 genes on the log TPM plot

locus.annotate <- function(z,a) {
  ifelse(rownames(z) == "CDKN2B-AS1", "CDKN2B-AS1",
         ifelse(rownames(z) == "MTAP", "MTAP",
                ifelse(rownames(z) == "DMRTA1", "DMRTA1",
                       ifelse(rownames(z) == "CDKN2A", "CDKN2A",
                              ifelse(rownames(z) == "CDKN2B", "CDKN2B",
                                     ifelse(rownames(z) == "MIR31HG", "MIR31HG",
                                            a))))))
}

#calculate logTPM for NT controls (all summed within replicates), and graph replicates

nt_bcs <- subset(barcode_key, guide_ctrl=="Scrambled control")
nt_bcs <- nt_bcs$bc

nt_genecounts1 <- subset(genecounts1, select = nt_bcs)
nt_genecounts2 <- subset(genecounts2, select = nt_bcs)
nt_genecounts <- as.data.frame(cbind(rowSums(nt_genecounts1), rowSums(nt_genecounts2)))

nt_genecounts_tpm <- nt_genecounts %>% mutate_each(logtpm)

nt_genecounts_tpm <- nt_genecounts_tpm %>%
  mutate(in_locus = locus.annotate(z=nt_genecounts_tpm, a=NA)
  )

nt_genecounts_tpm <- nt_genecounts_tpm %>%
  mutate(in_locus2 = locus.annotate(z=nt_genecounts_tpm, a=1)
  )

nt_ctrl_tpm <- ggplot(nt_genecounts_tpm %>% arrange(in_locus2), aes(x=V1, y=V2, color = in_locus, label = in_locus)) +
  geom_point() + geom_text_repel(vjust = 0, hjust = -0.8, show.legend=F, na.rm=T) +
  theme_bw() + xlab("log10(CPM+1)") + ylab("log10(CPM+1)") +
  labs(color = "9p21 Locus", title= "Non-targeting controls") + 
  scale_color_discrete(breaks=c("MTAP","CDKN2A","CDKN2B","CDKN2B-AS1","DMRTA1", "MIR31HG"))
nt_ctrl_tpm

ggsave(filename= file.path(save_path, "ctrl_nt_tpm.png"), plot=nt_ctrl_tpm, device=png, width = 6.5, height = 5)

#calculate logTPM for locus controls (all summed within replicates), and graph replicates

ctr_bcs <- subset(barcode_key, guide_ctrl=="Locus control")
ctr_bcs <- ctr_bcs$bc
#remove the missing barcode manually
#ctr_bcs <- ctr_bcs[-4]

ctr_genecounts1 <- subset(genecounts1, select = ctr_bcs)
ctr_genecounts2 <- subset(genecounts2, select = ctr_bcs)
ctr_genecounts <- as.data.frame(cbind(rowSums(ctr_genecounts1), rowSums(ctr_genecounts2)))

ctr_genecounts_tpm <- ctr_genecounts %>% mutate_each(logtpm)

ctr_genecounts_tpm <- ctr_genecounts_tpm %>%
  mutate(in_locus = locus.annotate(z=ctr_genecounts_tpm, a=NA)
  )

ctr_genecounts_tpm <- ctr_genecounts_tpm %>%
  mutate(in_locus2 = locus.annotate(z=ctr_genecounts_tpm, a=1)
  )

ctrl_tpm <- ggplot(ctr_genecounts_tpm %>% arrange(in_locus2), aes(x=V1, y=V2, color = in_locus, label = in_locus)) +
  geom_point() + geom_text_repel(vjust = 0, hjust = -1, show.legend=F, na.rm=T) +
  theme_bw() + xlab("log10(CPM+1)") + ylab("log10(CPM+1)") +
  labs(color = "9p21 Locus", title = "In-locus controls") +
  scale_color_discrete(breaks=c("MTAP","CDKN2A","CDKN2B","CDKN2B-AS1","DMRTA1", "MIR31HG"))
ctrl_tpm

ggsave(filename= file.path(save_path, "ctrl_locus_tpm.png"), plot=ctrl_tpm, device=png, width = 6.5, height = 5)


####extract TPM for 9p21 genes, plot TPM replicates together####

#get TPM values for each gene using the function assigned above

readcounts <- data.frame(read1 = colSums(genecounts1))
readcounts <- merge(readcounts, colSums(genecounts2), by="row.names", all=T)
colnames(readcounts) <- c("bc", "rep1", "rep2")

counts_tpm <- barcode_key

counts_rep1tpm <- barcode_key %>%
  mutate(rep_id = repnum1,
         seq_id = seqnum,
  )
counts_rep1tpm <- merge(counts_rep1tpm, colSums(genecounts1), by.x="bc",by.y="row.names", all=T)

counts_rep2tpm <- barcode_key %>%
  mutate(rep_id = repnum2,
         seq_id = seqnum,
  )
counts_rep2tpm <- merge(counts_rep2tpm, colSums(genecounts2), by.x="bc",by.y="row.names", all=T)

counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"CDKN2B-AS1",counts_tpm)
counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"MTAP",counts_tpm)
counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"CDKN2A",counts_tpm)
counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"CDKN2B",counts_tpm)
counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"DMRTA1",counts_tpm)
counts_tpm <- getcounts(genecounts1_tpm,genecounts2_tpm,"MIR31HG",counts_tpm)

counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "CDKN2B-AS1", counts_rep1tpm)
counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "MTAP", counts_rep1tpm)
counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "CDKN2A", counts_rep1tpm)
counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "CDKN2B", counts_rep1tpm)
counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "DMRTA1", counts_rep1tpm)
counts_rep1tpm <- getcounts_rep(genecounts1_tpm, "MIR31HG", counts_rep1tpm)

counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "CDKN2B-AS1", counts_rep2tpm)
counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "MTAP", counts_rep2tpm)
counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "CDKN2A", counts_rep2tpm)
counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "CDKN2B", counts_rep2tpm)
counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "DMRTA1", counts_rep2tpm)
counts_rep2tpm <- getcounts_rep(genecounts2_tpm, "MIR31HG", counts_rep2tpm)

write.csv(counts_tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_platemap_tpm.csv")
write.csv(counts_rep1tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_rep1_tpm.csv")
write.csv(counts_rep2tpm, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_rep2_tpm.csv")

#calculate linear regressions for each gene of interest
model_anril <-lm(`CDKN2B-AS1.y`~`CDKN2B-AS1.x`,data=counts_tpm)

model_mtap <-lm(MTAP.y~MTAP.x,data=counts_tpm)

model_cdkn2b <-lm(CDKN2B.y~CDKN2B.x,data=counts_tpm)

model_cdkn2a <-lm(CDKN2A.y~CDKN2A.x,data=counts_tpm)

model_dmrta1 <-lm(DMRTA1.y~DMRTA1.x,data=counts_tpm)

model_mir31hg <-lm(MIR31HG.y~MIR31HG.x,data=counts_tpm)

eqs <- list(
  eq_anril = lm_eq(model_anril),
  eq_mtap = lm_eq(model_mtap),
  eq_cdkn2b = lm_eq(model_cdkn2b),
  eq_cdkn2a = lm_eq(model_cdkn2a),
  eq_dmrta1 = lm_eq(model_dmrta1),
  eq_mir31hg = lm_eq(model_mir31hg)
)

#next line is just to check if lm_eq is working:
#eq_anril <- lm_eq(model_anril)

#plot CPM values included lm line

anril_reg <- ggplot(counts_tpm, aes(x=`CDKN2B-AS1.x`, y=`CDKN2B-AS1.y`)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("CDKN2B-AS1 CPM in 1st replicate") + ylab("CDKN2B-AS1 CPM in 2nd replicate")
anril_reg <- anril_reg + annotate(geom="text", x=55, y=50, label = as.call(eqs$eq_anril))
anril_reg

mtap_reg <- ggplot(counts_tpm, aes(x=MTAP.x, y=MTAP.y)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("MTAP CPM in 1st replicate") + ylab("MTAP CPM in 2nd replicate")
mtap_reg <- mtap_reg + annotate(geom="text", x=15, y=7, label = as.call(eqs$eq_mtap))
mtap_reg

cdkn2b_reg <- ggplot(counts_tpm, aes(x=CDKN2B.x, y=CDKN2B.y)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("CDKN2B CPM in 1st replicate") + ylab("CDKN2B CPM in 2nd replicate")
cdkn2b_reg <- cdkn2b_reg + annotate(geom="text", x=10, y=1, label = as.call(eqs$eq_cdkn2b))
cdkn2b_reg

cdkn2a_reg <- ggplot(counts_tpm, aes(x=CDKN2A.x, y=CDKN2A.y)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("CDKN2A CPM in 1st replicate") + ylab("CDKN2A CPM in 2nd replicate")
cdkn2a_reg <- cdkn2a_reg + annotate(geom="text", x=365, y=170, label = as.call(eqs$eq_cdkn2a))
cdkn2a_reg

dmrta1_reg <- ggplot(counts_tpm, aes(x=DMRTA1.x, y=DMRTA1.y)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("DMRTA1 CPM in 1st replicate") + ylab("DMRTA1 CPM in 2nd replicate")
dmrta1_reg <- dmrta1_reg + annotate(geom="text", x=4, y=0.3, label = as.call(eqs$eq_dmrta1))
dmrta1_reg

mir31hg_reg <- ggplot(counts_tpm, aes(x=MIR31HG.x, y=MIR31HG.y)) +
  geom_smooth(method='lm', se=F) + 
  geom_point(aes(color=guide_ctrl)) +
  #ylim(0,40)+xlim(0,40) +
  theme_bw() + xlab("MIR31HG CPM in 1st replicate") + ylab("MIR31HG CPM in 2nd replicate")
mir31hg_reg <- mir31hg_reg + annotate(geom="text", x=5, y=20, label = as.call(eqs$eq_dmrta1))
mir31hg_reg


ggsave(filename= file.path(save_path, "tpm_anril_reg.png"), plot=anril_reg, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "tpm_mtap_reg.png"), plot=mtap_reg, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "tpm_cdkn2a_reg.png"), plot=cdkn2a_reg, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "tpm_cdkn2b_reg.png"), plot=cdkn2b_reg, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "tpm_dmrta1_reg.png"), plot=dmrta1_reg, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "tpm_mir31hg_reg.png"), plot=mir31hg_reg, device=png, width = 6, height = 5)

####QC + compare read counts across replicates####

#extract read counts, %mito, and %ribo per well, and export as csv

readcounts <- data.frame(read1 = colSums(genecounts1))
readcounts <- merge(readcounts, colSums(genecounts2), by="row.names", all=T)
colnames(readcounts) <- c("bc", "rep1", "rep2")

mt1 <- data.frame(genecounts1[grep('^MT-',rownames(genecounts1)),])
mt2 <- data.frame(genecounts2[grep('^MT-',rownames(genecounts2)),])
mt_counts <- data.frame(mt1 = colSums(mt1))
mt_counts <- merge(mt_counts, colSums(mt2), all=T, by="row.names")
colnames(mt_counts) <- c("bc", "mt_rep1", "mt_rep2")

ribo1 <- data.frame(genecounts1[grep('^RPS|^RPL',rownames(genecounts1)),])
ribo2 <- data.frame(genecounts2[grep('^RPS|^RPL',rownames(genecounts2)),])
ribo_counts <- data.frame(ribo1 = colSums(ribo1))
ribo_counts <- merge(ribo_counts, colSums(ribo2), all=T, by="row.names")
colnames(ribo_counts) <- c("bc", "ribo_rep1", "ribo_rep2")


readcounts <- merge(readcounts, mt_counts, by="bc")
readcounts <- merge(readcounts, ribo_counts, by="bc")
readcounts <- readcounts %>% mutate(
  mtper1 = rep1/mt_rep1,
  mtper2 = rep2/mt_rep2,
  riboper1 = rep1/ribo_rep1,
  riboper2 = rep2/ribo_rep2
)
readcounts <- merge(barcode_key, readcounts, by="bc", all=T)

#write.csv(readcounts, "/Volumes/broad_mcl/members_dir/bschmand/data/mac_seq/cfb_dxk99/dxk99_cfb_qc_counts.csv")


#plot a bunch of QC metrics
#plot: read counts by rep 1 and 2

#calculate regression line for read counts (added for NNFC 2024)
model_counts <-lm(rep1~rep2,data=readcounts)
counts_eq <- lm_eq(model_counts)

count_scatter <- ggplot(data=readcounts, aes(x=rep1, y=rep2)) +
  geom_point(aes(color=guide_ctrl)) +
  geom_smooth(method='lm', se=F) +
  theme_bw() + xlab("Replicate 1 read counts") + ylab("Replicate 2 read counts")
count_scatter <- count_scatter + annotate(geom="text", x=2000000, y=750000, label = as.call(counts_eq))
count_scatter

#plot: mitochondrial read counts by replicate
mt_scatter <- ggplot(readcounts, aes(x=mt_rep1, y=mt_rep2, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 1 mitochondrial reads") + ylab("Replicate 2 mitochondrial reads")
mt_scatter

ribo_scatter <- ggplot(readcounts, aes(x=ribo_rep1, y=ribo_rep2, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 1 ribosomal reads") + ylab("Replicate 2 ribosomal reads")
ribo_scatter

#plot: mito fraction by replicate 
mt_pscatter <- ggplot(readcounts, aes(x=mtper1, y=mtper2, color=guide_ctrl)) +
  geom_point() +
  ylim(0,35) + xlim(0,35) +
  theme_bw() + xlab("Replicate 1 %mt") + ylab("Replicate 2 %mt")
mt_pscatter

#for 1 replicate: plot read depth by %mito (see how mito% changes by depth)
mt_pscatter <- ggplot(readcounts, aes(x=rep1, y=mtper1, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 1 read count") + ylab("Replicate 1 %mt")
mt_pscatter

ggsave(filename= file.path(save_path, "qc_mt_pct_scatter1.png"), plot=mt_pscatter, device=png, width = 6, height = 5)

#for 1 replicate: plot read depth by %mito (see how mito% changes by depth)
ribo_pscatter <- ggplot(readcounts, aes(x=rep1, y=riboper1, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 1 read count") + ylab("Replicate 1 %ribo")
ribo_pscatter

ggsave(filename= file.path(save_path, "qc_mt_pct_scatter1.png"), plot=mt_pscatter, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "qc_ribo_pct_scatter1.png"), plot=ribo_pscatter, device=png, width = 6, height = 5)

#for 2nd replicate: plot read depth by %mito (see how mito% changes by depth)
mt_pscatter <- ggplot(readcounts, aes(x=rep2, y=mtper2, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 2 read count") + ylab("Replicate 2 %mt")
mt_pscatter

ribo_pscatter <- ggplot(readcounts, aes(x=rep2, y=riboper2, color=guide_ctrl)) +
  geom_point() +
  theme_bw() + xlab("Replicate 2 read count") + ylab("Replicate 2 %ribo")
ribo_pscatter

#plot TSS effectiveness

ggsave(filename= file.path(save_path, "qc_count_scatter.png"), plot=count_scatter, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "qc_mt_scatter.png"), plot=mt_scatter, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "qc_ribo_scatter.png"), plot=ribo_scatter, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "qc_mt_pct_scatter2.png"), plot=mt_pscatter, device=png, width = 6, height = 5)
ggsave(filename= file.path(save_path, "qc_ribo_pct_scatter2.png"), plot=ribo_pscatter, device=png, width = 6, height = 5)






