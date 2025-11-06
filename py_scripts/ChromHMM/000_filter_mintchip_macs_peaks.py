import csv
import json
import math
import os
import sys
import pandas as pd
import statistics
import numpy as np

"""
Cleaning macs3 peaks based on fold change and width
"""


def create_dir(dir_path):
    if not os.path.exists(dir_path):
        os.system("mkdir %s" % dir_path)


params = json.load(open("ChromHMM_params.json", "r"))
report = csv.writer(open(os.path.join("../report/macs_stats_minchip.csv"), "w"))
report.writerow(["CellType", "Stimulation", "histone_mark", "Mean fc", "SD fc", "25th fc", "75th fc", "Num initial",
                 "Num filtered"])

mint_cip_root = params["root_mintchip_data"]
for tissue in os.listdir(mint_cip_root):
    tissue_path = os.path.join(mint_cip_root, tissue)
    for stimulation in os.listdir(tissue_path):
        if stimulation == ".DS_Store" or stimulation != "basal":
            continue
        stimulation_path = os.path.join(tissue_path, stimulation)
        for histone_mark in os.listdir(stimulation_path):
            if histone_mark == ".DS_Store" or histone_mark == "WCE":
                continue
            histone_path = os.path.join(stimulation_path, histone_mark)
            file_name = "%s_WCE" % histone_mark
            peaks_dir_input = os.path.join(histone_path, "merged_macs", file_name)
            if not os.path.exists(os.path.join(peaks_dir_input, "%s_narrow_peaks.narrowPeak" % file_name)):
                print("Missing peak files; %s" % os.path.join(peaks_dir_input, "%s_narrow_peaks.narrowPeak" % file_name))
                continue
            input_data = pd.read_csv(os.path.join(peaks_dir_input, "%s_narrow_peaks.narrowPeak" % file_name),
                                     sep='\t', header=None)
            output_path = os.path.join(peaks_dir_input, "%s_narrow_peaks_filtered.narrowPeak" % file_name)

            fold_change = input_data[6]
            len_init = len(input_data)

            iqr25_fc = np.percentile(fold_change, 25, interpolation='midpoint')
            iqr75_width = np.percentile(input_data[4], 75, interpolation='midpoint')
            median = statistics.median(input_data[6])
            # input_data = input_data[input_data[6] > iqr25]
            input_data = input_data[input_data[6] > median]
            input_data = input_data[input_data[4] > iqr75_width]
            input_data = input_data[input_data[7] > -math.log10(.05)]

            input_data.to_csv(output_path, sep='\t', header=None, index=False)
            report.writerow([tissue, stimulation, histone_mark,
                                                           statistics.mean(fold_change),
                                                           statistics.stdev(fold_change),
                                                           np.percentile(fold_change, 25, interpolation='midpoint'),
                                                           np.percentile(fold_change, 75, interpolation='midpoint'),
                                                           len_init, len(input_data)])
