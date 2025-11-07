import json
import os


def get_uger_path(path):
    return path.replace("/broad/", "/broad/")


def create_dir(path):
    if not os.path.exists(path):
        os.system("mkdir %s" % path)


params = json.load(open("ChromHMM_params.json", "r"))
root_binaries = "9p21/ChromHMM/binerized/jointly_bed/"
root_models = "9p21/ChromHMM/models/"
create_dir(root_models)

qsub_individual_path = os.path.join(params["qsub_individual_path"], "chromhmm_modeling")
create_dir(qsub_individual_path)

scripts_paths = []


def iterate_on_num_states():
    joint_model_path = root_binaries
    for num_states in [10, 12, 14, 15, 16, 18, 20, 22]:
        output_model_dir = os.path.join(root_models, "%d" % num_states)
        create_dir(output_model_dir)

        commands = ["INPUT_DIR=%s" % get_uger_path(joint_model_path),
                    "OUTPUT_DIR=%s/" % get_uger_path(output_model_dir),
                    "NUM_STATS=%d" % num_states,
                    "ASSEMBLY=hg38", "", "",
                    "java -mx16000M -jar $ChromHMM LearnModel $INPUT_DIR $OUTPUT_DIR $NUM_STATS $ASSEMBLY",
                    ""]
        out_file_name = "%s_model_%d.sh" % ("chromhmm_modeling", num_states)
        out_file_path = os.path.join(qsub_individual_path, out_file_name)
        scripts_paths.append(get_uger_path(out_file_path))
        fout = open(out_file_path, "w")
        fout.write(params["qsub_header"])
        fout.write("\n".join(commands))
        fout.close()

iterate_on_num_states()
print("%d py_r_scripts were generated" % len(scripts_paths))
qsub_output_path = os.path.join(params["qsub_output_path"], "chromhmm_modeling")
create_dir(qsub_output_path)

call_sh_path = os.path.join(params["qsub_call_path"], "qsub_chromhmm_modeling.sh")
fout = open(call_sh_path, "w")
fout.write(params["qsub_header"])
for a_script in scripts_paths:
    fout.write(
        "qsub -l h_vmem=32G -l h_rt=6:00:00 -o %s -e %s %s\n\n" % (get_uger_path(qsub_output_path),
                                                                   get_uger_path(qsub_output_path), a_script))
fout.write("qstat\n")
fout.close()

os.system("chmod 777 %s" % call_sh_path)