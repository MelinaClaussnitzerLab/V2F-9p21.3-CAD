import os
import json


def get_uger_path(path):
    return path.replace("/broad/", "/broad/")


def create_dir(path):
    if not os.path.exists(path):
        os.system("mkdir %s" % path)


def get_cell_type_abbreviation(tissue):
    if tissue in ["Cardiac_fibroblasts_atrial", "Cardiac_Fibroblasts_atrial_May_11_2021"]:
        return "CFA"
    if tissue in ["Endothelial_coronary_artery", "Endothelial_coronary_artery_May_11_2021"]:
        return "ECA"
    if tissue in ["VSMC_coronary_artery", "VSMC_coronary_artery_May_4_2022"]:
        return "VSMCCA"
    if tissue in ["pericytes_cortical", "Pericytes_cortical_May_4_2022"]:
        return "PC"
    if tissue in ["VSMC_pulmonary_artery", "VSMC_pulmonary_artery_May_11_2021"]:
        return "VSMCPA"
    if tissue in ["Pericytes_adipose", "Pericytes_adipose_May_11_2021"]:
        return "PA"

    if tissue in ["Pericytes_adipose_May_4_2022", "Cardiac_fibroblasts_ventricular"]:
        return "NA"
    return "NA"


params = json.load(open("ChromHMM_params.json", "r"))
atac_root = params["root_atac_data"]
chip_root = params["root_mintchip_data"]
root_output = "9p21/ChromHMM/binerized/jointly_bed/"

cellmarkfiletable = open(os.path.join(root_output, "cellmarkfiletable.txt"), "w")

for tissue in os.listdir(atac_root):
    tissue_abb = get_cell_type_abbreviation(tissue)
    if tissue == ".DS_Store" or tissue_abb == "NA":
        continue
    tissue_path = os.path.join(atac_root, tissue)
    for stimulation in os.listdir(tissue_path):
        if stimulation == ".DS_Store" or stimulation != "basal":  # this is because we use basal for ChromHMM
            continue
        stimulation_path = os.path.join(tissue_path, stimulation, "macs3_peaks_hg38")
        for a_file in os.listdir(stimulation_path):
            if "_filtered.narrowPeak" not in a_file:
                continue
            source = os.path.join(stimulation_path, a_file)
            output_name = "%s_%s" % (tissue_abb, a_file.replace("_filtered.narrowPeak", ".bed"))
            destination = os.path.join(root_output, output_name)
            os.system("cp %s %s" % (source, destination))
            cellmarkfiletable.write("%s\tATAC\t%s\n" % (tissue_abb, output_name))

for tissue in os.listdir(chip_root):
    tissue_abb = get_cell_type_abbreviation(tissue)
    if tissue == ".DS_Store" or tissue_abb == "NA":
        continue
    tissue_path = os.path.join(chip_root, tissue)
    for stimulation in os.listdir(tissue_path):
        if stimulation == ".DS_Store" or stimulation != "basal":
            continue
        stimulation_path = os.path.join(tissue_path, stimulation)
        for histone_mark in os.listdir(stimulation_path):
            if histone_mark == ".DS_Store" or histone_mark == "WCE":
                continue
            histone_path = os.path.join(stimulation_path, histone_mark)
            file_name = "%s_WCE" % histone_mark
            source = os.path.join(histone_path, "merged_macs", file_name,
                                  "%s_narrow_peaks_filtered.narrowPeak" % file_name)
            output_name = "%s_%s.bed" % (tissue_abb, histone_mark)
            destination = os.path.join(root_output, output_name)
            os.system("cp %s %s" % (source, destination))
            cellmarkfiletable.write("%s\t%s\t%s\n" % (tissue_abb, histone_mark, output_name))


cellmarkfiletable.close()

commands = ["DATA_dir=%s" % get_uger_path(root_output), "java -mx16000M -jar $ChromHMM BinarizeBed -b 200 -peaks \\",
            "%s \\" % params["chromsize_hg38"], "$DATA_dir \\", "$DATA_dir/cellmarkfiletable.txt \\", "$DATA_dir", "",
            ""]

scripts_paths = []
qsub_individual_path = os.path.join(params["qsub_individual_path"], "chromhmm_binerize")
create_dir(qsub_individual_path)
out_file_name = "chromhmm_joint_modeling_binerize_bed.sh"
out_file_path = os.path.join(qsub_individual_path, out_file_name)
scripts_paths.append(get_uger_path(out_file_path))
fout = open(out_file_path, "w")
fout.write(params["qsub_header"])
fout.write("\n".join(commands))
fout.close()

print("%d scripts were generated" % len(scripts_paths))
qsub_output_path = os.path.join(params["qsub_output_path"], "chromhmm_binerize")
create_dir(qsub_output_path)

call_sh_path = os.path.join(params["qsub_call_path"], "qsub_%s.sh" % "chromhmm_binerize")
fout = open(call_sh_path, "w")
fout.write(params["qsub_header"])
for a_script in scripts_paths:
    fout.write(
        "qsub -l h_vmem=16G -l h_rt=1:00:00 -o %s -e %s %s\n\n" % (get_uger_path(qsub_output_path),
                                                                   get_uger_path(qsub_output_path), a_script))
fout.write("qstat\n")
fout.close()

os.system("chmod 777 %s" % call_sh_path)