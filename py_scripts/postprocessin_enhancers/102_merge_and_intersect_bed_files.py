import os
import pandas as pd


roadmap_enhancers = pd.read_csv("../data/roadmap_enhancers.csv")
chromhmm_enhancers = pd.read_csv("../data/ChromHMM_JointModel_enhances.csv")

print(roadmap_enhancers.columns)
print(chromhmm_enhancers.columns)
roadmap_enhancers["NAME"] = roadmap_enhancers["NAME"]+";"+roadmap_enhancers["fname"]
chromhmm_enhancers["Name"] = chromhmm_enhancers["Name"]+";"+chromhmm_enhancers["Fname"]

roadmap_enhancers = roadmap_enhancers.drop(columns=["fname"])
chromhmm_enhancers = chromhmm_enhancers.drop(columns=["Fname"])

roadmap_enhancers.to_csv("../data/roadmap_enhancers.bed", sep="\t", header=None, index=False)
chromhmm_enhancers.to_csv("../data/chromHMM_enhancers.bed", sep="\t", header=None, index=False)


command = ["use Anaconda3", "ROOT=conda_envs/rose_env",
           "source activate $ROOT", "", ""]

root = "postProcessing_enhancers/data/"
intersect_output = os.path.join(root, "ChromHMM_JointModel_roadmap_intersect.bed")
roadmap_enhancers_orig_path = os.path.join(root, "roadmap_enhancers.bed")
roadmap_enhancers_sorted_path = os.path.join(root, "roadmap_enhancers_sorted.bed")
roadmap_enhancers_merged_path = os.path.join(root, "roadmap_enhancers_merged.bed")


chromhmm_enhancers_orig_path = os.path.join(root, "chromHMM_enhancers.bed")
chromhmm_enhancers_sorted_path = os.path.join(root, "chromHMM_enhancers_sorted.bed")
chromhmm_enhancers_merged_path = os.path.join(root, "chromHMM_enhancers_merged.bed")


command.append("sort -k1,1 -k2,2n %s > %s" % (roadmap_enhancers_orig_path, roadmap_enhancers_sorted_path))
command.append("sort -k1,1 -k2,2n %s > %s" % (chromhmm_enhancers_orig_path, chromhmm_enhancers_sorted_path))
command.append("bedtools merge -i %s > %s" % (roadmap_enhancers_sorted_path, roadmap_enhancers_merged_path))
command.append("bedtools merge -i %s > %s" % (chromhmm_enhancers_sorted_path, chromhmm_enhancers_merged_path))
command.append("bedtools intersect -wa -wb -a %s -b %s -sorted > %s " %
               (chromhmm_enhancers_merged_path, roadmap_enhancers_merged_path, intersect_output))

fout = open("../scripts/intersect_chromHMM_enhancer.sh", "w")
fout.write("\n".join(command))
fout.close()

