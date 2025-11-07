import csv

import pandas as pd
import os


output = "../data/ChromHMM_JointModel_enhances.csv"
root = "/Volumes/broad_mcl/database/9p21/ChromHMM/run/models/18"
enhancer_annotations = ["E5", "E6", "E7", "E8", "E10"]
file_names = ["VSMCCA_18_segments.bed", "VSMCPA_18_segments.bed", "CFA_18_segments.bed",
              "PA_18_segments.bed", "ECA_18_segments.bed"]
writer = csv.writer(open(output, "w"))
writer.writerow(["Fname", "CHR", "START", "END", "Name"])
for a_file in file_names:
    data = pd.read_csv(os.path.join(root, a_file), sep="\t", header=None, names=["CHR", "START", "END", "Name"])
    for an_enhancer in enhancer_annotations:
        curr_data = data[data["Name"] == an_enhancer]
        for index, row in curr_data.iterrows():
            out_row = [a_file, row["CHR"], row["START"], row["END"], row["Name"]]
            writer.writerow(out_row)