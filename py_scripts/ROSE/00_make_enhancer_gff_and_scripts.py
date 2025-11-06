import os

def create_dir(path):
    if not os.path.exists(path):
        os.system("mkdir %s" % path)


root_call_qsub = "qsub_files/qsub_call_scripts/"
root_output_qsub = "qsub_files/qsub_outputs/rose_9p21/"
root_indi_qsub = "qsub_files/qsub_individual_scripts/rose_9p21/"
create_dir(root_output_qsub)
create_dir(root_indi_qsub)

root_chromhmm_results = "chromhmm/models/w_atac_joint_merged_mintchip_beds/15_annotated"
enhancer_annotations = ["07_EnhAc", "08_EnhAc", "09_Enh", "10_EnhP"]
region_of_interest = {"chr": "chr9", "start": 21171056, "end": 22400575}

root_output = "ROSE/Run_on_hg19"
bam_paths = {
    "VSMC":
        {
            "RANKING_BAM": "VSMC_pulmonary_artery_May_11_2021/basal/H3K27ac/bam_merged_to_one/H3K27ac.sorted.bam",
            "CONTROL_BAM": "VSMC_pulmonary_artery_May_11_2021/basal/WCE/bam_merged_to_one/WCE.sorted.bam"
        },
    "CFA":
        {
            "RANKING_BAM": "Cardiac_Fibroblasts_atrial_May_11_2021/basal/H3K27ac/bam_merged_to_one/H3K27ac.sorted.bam",
            "CONTROL_BAM": "Cardiac_Fibroblasts_atrial_May_11_2021/basal/WCE/bam_merged_to_one/WCE.sorted.bam"
        },
    "PA":
        {
            "RANKING_BAM": "Pericytes_adipose_May_11_2021/basal/H3K27ac/bam_merged_to_one/H3K27ac.sorted.bam",
            "CONTROL_BAM": "Pericytes_adipose_May_11_2021/basal/WCE/bam_merged_to_one/WCE.sorted.bam"
        },
    "ECA":
        {
            "RANKING_BAM": "Endothelial_coronary_artery_May_11_2021/basal/H3K27ac/bam_merged_to_one/H3K27ac.sorted.bam",
            "CONTROL_BAM": "Endothelial_coronary_artery_May_11_2021/basal/WCE/bam_merged_to_one/WCE.sorted.bam"
        }
}
scripts_list = []
for a_file in ["VSMC.bed", "CFA.bed", "PA.bed", "ECA.bed"]:
    for location in ["genome_wide", "9p21"]:
        tissue_name = a_file.replace(".bed", "")
        curr_output = os.path.join(root_output, tissue_name)
        create_dir(curr_output)
        curr_output = os.path.join(curr_output, location)
        create_dir(curr_output)
        result_path = os.path.join(curr_output, "results")
        create_dir(result_path)
        gff_file_path = os.path.join(curr_output, "%s_%s.gff" % (tissue_name, location))
        fout = open(gff_file_path, "w")
        fin = open(os.path.join(root_chromhmm_results, a_file), "r")
        counter = 0
        for a_line in fin:
            for an_enhancer in enhancer_annotations:
                if an_enhancer in a_line:
                    # chr10	119200	119600	03_Prom	0	.	119200	119600	239,57,31
                    content = a_line.split("\t")
                    chr = content[0]
                    if len(chr) > 5:
                        continue
                    start = content[1]
                    end = content[2]
                    if location == "9p21":
                        if chr != region_of_interest["chr"]:
                            continue
                        if not (region_of_interest["start"] <= int(start) <= region_of_interest["end"] or
                                region_of_interest["start"] <= int(end) <= region_of_interest["end"]):
                            continue

                    """
                    1: chromosome (chr#)
                    2: unique ID for each constituent enhancer region
                    4: start of constituent
                    5: end of constituent
                    7: strand (+,-,.)
                    9: unique ID for each constituent enhancer region
                    """
                    counter += 1
                    name = "l%s_%d" % (an_enhancer, counter)
                    ID = "%d" % counter
                    fout.write("%s	%s		%s	%s		.		%s\n"
                               % (chr, name, start, end, name))
        fout.close()
        fin.close()
        script_name = "%s_%s.sh" % (tissue_name, location)
        script_path = os.path.join(root_indi_qsub, script_name)
        fout = open(script_path, "w")
        fout.write("#!/bin/bash\nsource /broad/software/scripts/useuse\n\n")
        fout.write("use Anaconda3\nsource activate conda_envs/rose_env\n\n")
        fout.write("PATHTO=rose/scripts/\n\n")
        fout.write("cd $PATHTO\n")
        fout.write("PYTHONPATH=$PATHTO/lib\nexport PYTHONPATH\nexport PATH=$PATH:$PATHTO/bin\n\n")
        fout.write("INPUT_CONSTITUENT_GFF=%s\n" % gff_file_path)
        fout.write("RANKING_BAM=%s\n" % bam_paths[tissue_name]["RANKING_BAM"])
        fout.write("CONTROL_BAM=%s\n" % bam_paths[tissue_name]["CONTROL_BAM"])
        fout.write("OUTPUT_DIRECTORY=%s\n\n" % result_path)
        fout.write("ROSE_main.py -g HG19 -i $INPUT_CONSTITUENT_GFF -r $RANKING_BAM -c $CONTROL_BAM -o $OUTPUT_DIRECTORY")
        fout.write("\n\n")
        scripts_list.append(script_path)


call_script_name = os.path.join(root_call_qsub, "rose_9p21_hg19.sh")
fout = open(call_script_name, "w")
counter = 0
qsub_memory = 32
qsub_hours = 12
for _ in scripts_list:
    counter += 1
    fout.write('echo "%d"\n' % counter)
    fout.write("qsub -l h_vmem=%dG -l h_rt=%d:00:00 -o %s -e %s %s\n\n" % (
        qsub_memory, qsub_hours, root_output_qsub, root_output_qsub, _))
fout.write("qstat\n\n")
fout.close()
os.system("chmod 777 %s" % call_script_name)
print(["# scripts:", len(scripts_list)])