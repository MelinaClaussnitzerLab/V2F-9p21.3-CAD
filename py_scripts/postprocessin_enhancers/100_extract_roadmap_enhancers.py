import csv
import os
import pandas as pd


root = "../data/all.mnemonics.bedFiles"
writer = csv.writer(open("../data/roadmap_enhancers.csv", "w"))
writer.writerow(["fname", "CHR", "START", "END", "NAME"])
for a_file in os.listdir(root):
    if a_file == ".DS_Store":
        continue
    print(a_file)
    data = pd.read_csv(os.path.join(root, a_file), compression="gzip", sep="\t", header=None,
                       names=["CHR", "START", "END", "NAME"])
    for an_enhancer in ["7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk"]:
        curr_data = data[data["NAME"] == an_enhancer]
        for index, row in curr_data.iterrows():
            writer.writerow([a_file, row["CHR"], row["START"], row["END"], row["NAME"]])

    """
    7_EnhG1	Genic enhancer1	GreenYellow	194,225,5
    8_EnhG2	Genic enhancer2	GreenYellow	194,225,5
    9_EnhA1	Active Enhancer 1	Orange	255,195,77
    10_EnhA2	Active Enhancer 2	Orange	255,195,77
    11_EnhWk	Weak Enhancer	Yellow	255,255,0
    """

