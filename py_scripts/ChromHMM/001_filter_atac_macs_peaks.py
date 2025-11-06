import json
import os
import csv
import pandas as pd
import statistics
import numpy as np
import math


def create_dir(dir_path):
    if not os.path.exists(dir_path):
        os.system("mkdir %s" % dir_path)


report = csv.writer(open(os.path.join("../report/macs_stats_atac.csv"), "w"))
report.writerow(["CellType", "Stimulation", "histone_mark", "Mean fc", "SD fc", "25th fc", "75th fc", "Num initial",
                 "Num filtered"])
params = json.load(open("ChromHMM_params.json", "r"))
data_folder = params["root_atac_data"]
for tissue in os.listdir(data_folder):
    if tissue == ".DS_Store":
        continue
    tissue_path = os.path.join(data_folder, tissue)
    for stimulation in os.listdir(tissue_path):
        if stimulation == ".DS_Store" or stimulation != "basal":  # this is because we use basal for ChromHMM
            continue
        stimulation_path = os.path.join(tissue_path, stimulation)
        macs3_path = os.path.join(stimulation_path, "macs3_peaks_hg38")

        def process_a_file(input_data):
            fold_change = input_data[6]
            len_init = len(input_data)

            iqr25 = np.percentile(fold_change, 25, interpolation='midpoint')
            iqr75_width = np.percentile(input_data[4], 75, interpolation='midpoint')
            median = statistics.median(input_data[6])

            input_data = input_data[input_data[6] > median]
            input_data = input_data[input_data[4] > iqr75_width]
            input_data = input_data[input_data[7] > -math.log10(.05)]

            report.writerow([tissue, stimulation, "ATAC",
                             statistics.mean(fold_change),
                             statistics.stdev(fold_change),
                             np.percentile(fold_change, 25, interpolation='midpoint'),
                             np.percentile(fold_change, 75, interpolation='midpoint'),
                             len_init, len(input_data)])
            return input_data


        basal_peak_in_path = os.path.join(macs3_path, "basal_peaks.narrowPeak")
        basal_peak_out_path = os.path.join(macs3_path, "basal_peaks_filtered.narrowPeak")
        basal_peak = pd.read_csv(basal_peak_in_path,
                                 sep='\t', header=None)

        input_data = process_a_file(basal_peak)
        input_data.to_csv(basal_peak_out_path, sep='\t', header=None, index=False)

        basal_peak_rep_in_path = os.path.join(macs3_path, "basal_replicate_peaks.narrowPeak")
        if os.path.exists(basal_peak_rep_in_path):
            basal_peak_rep_out_path = os.path.join(macs3_path, "basal_replicate_peaks_filtered.narrowPeak")
            basal_peak_rep = pd.read_csv(basal_peak_rep_in_path,
                                         sep='\t', header=None)

            input_data = process_a_file(basal_peak_rep)
            input_data.to_csv(basal_peak_rep_out_path, sep='\t', header=None, index=False)
