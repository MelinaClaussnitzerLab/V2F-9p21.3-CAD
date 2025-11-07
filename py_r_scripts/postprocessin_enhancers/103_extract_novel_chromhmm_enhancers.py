import csv
import os

import pandas as pd

root = "postProcessing_enhancers/data/"
writer = csv.writer(open(os.path.join(root, "novel_chromhmm_enhancers.bed"), "w"), delimiter="\t")
chrom_hmm_peaks = pd.read_csv(os.path.join(root, "chromHMM_enhancers_merged.bed"), sep="\t",
                              header=None, names=["CHR", "START", "END"])
roadmap_enhancers_peaks = pd.read_csv(os.path.join(root, "roadmap_enhancers_merged.bed"), sep="\t",
                              header=None, names=["CHR", "START", "END"])

matched_paeks = pd.read_csv(os.path.join(root, "ChromHMM_JointModel_roadmap_intersect.bed"), sep="\t",
                            header=None, names=["CHR", "START", "END", "chr2", "start2", "end2"])

Num_chrom_hmm_peaks = len(chrom_hmm_peaks)
Num_roadmap_enhancer_peaks = len(roadmap_enhancers_peaks)
Num_chromhmm_roadmap_intersecting = 0
Num_chromhmm_unique = 0
for index, row in chrom_hmm_peaks.iterrows():
    chr = row["CHR"]
    start = row["START"]
    end = row["END"]
    curr_matched_paeks = matched_paeks[matched_paeks["CHR"] == chr]
    curr_matched_paeks = curr_matched_paeks[curr_matched_paeks["START"] == start]
    curr_matched_paeks = curr_matched_paeks[curr_matched_paeks["END"] == end]
    if len(curr_matched_paeks) == 0:
        writer.writerow([chr, start, end])
        Num_chromhmm_unique += 1
    else:
        Num_chromhmm_roadmap_intersecting += 1

Num_roadmap_chromhmm_intersecting = 0
Num_roadmap_unique = 0
for index, row in roadmap_enhancers_peaks.iterrows():
    chr = row["CHR"]
    start = row["START"]
    end = row["END"]
    curr_matched_paeks = matched_paeks[matched_paeks["chr2"] == chr]
    curr_matched_paeks = curr_matched_paeks[curr_matched_paeks["start2"] == start]
    curr_matched_paeks = curr_matched_paeks[curr_matched_paeks["end2"] == end]
    if len(curr_matched_paeks) == 0:
        Num_roadmap_unique += 1
    else:
        Num_roadmap_chromhmm_intersecting += 1

writer = csv.writer(open(os.path.join(root, "chromhmm_vs_roadmap_counts_Feb20_2025.csv"), "w"))
writer.writerow(["Num_chrom_hmm_peaks", "Num_roadmap_enhancer_peaks",
                 "Num_chromhmm_roadmap_intersecting", "Num_chromhmm_unique",
                 "Num_roadmap_chromhmm_intersecting", "Num_roadmap_unique"])
writer.writerow([Num_chrom_hmm_peaks, Num_roadmap_enhancer_peaks,
                 Num_chromhmm_roadmap_intersecting, Num_chromhmm_unique,
                 Num_roadmap_chromhmm_intersecting, Num_roadmap_unique])
