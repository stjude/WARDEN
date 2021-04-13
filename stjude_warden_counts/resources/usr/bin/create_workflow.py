#!/usr/bin/env python
# warden_python 1.0.0


import sys
import dxpy
import re
import json

parameter_file = sys.argv[1]
count_samples_file = sys.argv[2]
final_sample_list_id = sys.argv[3]
comparisons_limma_id = sys.argv[4]
comparisons_all_id = sys.argv[5]

APPS = open("/app_data/script_names.txt")
app_names = {}
for line in APPS:
    line = line.rstrip("\n")
    var, val = line.split("\t")
    app_names[var] = val

# parse parameter file
parameters = {}
parameters_list = ""
PAR = open(parameter_file)
for line in PAR:
    parameters_list += line
    line = line.rstrip("\n")
    var, val = line.split("\t")
    parameters[var] = val

samples = {}
COUNTFILES = open(count_samples_file)
for line in COUNTFILES:
    line = line.rstrip("\n")
    sample_name, countfile = line.split("\t")
    samples[sample_name] = countfile


def build_workflow():
    if parameters["folder_provided"] == "false":
        wf = dxpy.new_dxworkflow(name='WARDEN_workflow',
                                 description='RNA-SEQ Workflow',
                                 output_folder=parameters["Output"],
                                 )
    else:
        wf = dxpy.new_dxworkflow(name='WARDEN_workflow',
                                 description='RNA-SEQ Workflow',
                                 )
    wf_outputs = []

    combine_counts_applet = dxpy.search.find_one_data_object(classname="applet",
                                                             name=app_names["combine_counts"],
                                                             state="closed",
                                                             return_handler=True)
    limma_applet = dxpy.search.find_one_data_object(classname="applet",
                                                    name=app_names["limma"],
                                                    state="closed",
                                                    return_handler=True)
    simple_DE_applet = dxpy.search.find_one_data_object(classname="applet",
                                                        name=app_names["simple_DE"],
                                                        state="closed",
                                                        return_handler=True)

    sample_num = 0
    htseq_results = [dxpy.dxlink(count_id) for count_id in samples.values()]

    combine_input = {"count_files": htseq_results, "name_value": "htseq", "sample_files": [dxpy.dxlink(final_sample_list_id)]}
    combine_counts_stage_id = wf.add_stage(combine_counts_applet, stage_input=combine_input, instance_type="azure:mem2_ssd1_x1", name="COMBINE HTSEQ")
    wf_outputs += [
        {
            "name": "combined_counts",
            "class": "file",
            "outputSource": {
                "$dnanexus_link": {
                    "stage": combine_counts_stage_id,
                    "outputField": "count_file"
                }
            }
        },
    ]

    if parameters["limma_DE_viewer"] != "None":
        limma_viewer_project, limma_viewer_file = parameters["limma_DE_viewer"].split(":")
        limma_viewer_link = dxpy.dxlink({"project": limma_viewer_project, "id": limma_viewer_file})

    if parameters["limma_runnable"] == "true":
        limma_input = {"input_count_file": dxpy.dxlink({"stage": combine_counts_stage_id, "outputField": "count_file"}), "sample_list_file": dxpy.dxlink(final_sample_list_id), "calcNormFactors_method": parameters["calcNormFactors_method"], "filter_count_type": parameters["filter_count_type"], "filter_count": int(parameters["filter_count"]), "p_value_adjust": parameters["p_value_adjust"], "contrasts_file": dxpy.dxlink(comparisons_limma_id)}
        if parameters["limma_DE_viewer"] != "None":
            limma_input["difex_viewer"] = limma_viewer_link
        limma_stage_id = wf.add_stage(limma_applet, stage_input=limma_input, instance_type="azure:mem1_ssd1_x4", name="LIMMA")
        wf_outputs += [
            {
                "name": "limma_outfiles",
                "class": "array:file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": limma_stage_id,
                        "outputField": "out_files"
                    }
                }
            },
            {
                "name": "limma_viewer",
                "class": "record",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": limma_stage_id,
                        "outputField": "viewer_bookmark"
                    }
                }
            },
        ]
    simple_DE_input = {"input_count_file": dxpy.dxlink({"stage": combine_counts_stage_id, "outputField": "count_file"}), "sample_list_file": dxpy.dxlink(final_sample_list_id), "contrasts_file": dxpy.dxlink(comparisons_all_id)}
    if parameters["limma_DE_viewer"] != "None":
        simple_DE_input["difex_viewer"] = limma_viewer_link
    simple_DE_stage_id = wf.add_stage(simple_DE_applet, stage_input=simple_DE_input, instance_type="azure:mem1_ssd1_x4", name="SIMPLE DIFFERENTIAL_EXPRESSION")
    wf_outputs += [
        {
            "name": "simple_DE_outfiles",
            "class": "array:file",
            "outputSource": {
                "$dnanexus_link": {
                    "stage": simple_DE_stage_id,
                    "outputField": "out_files"
                }
            }
        },
        {
            "name": "simple_DE_viewer",
            "class": "record",
            "outputSource": {
                "$dnanexus_link": {
                    "stage": simple_DE_stage_id,
                    "outputField": "viewer_bookmark"
                }
            }
        },
    ]

    wf.update(workflow_outputs=wf_outputs)
    wf.close()
    return wf.get_id()


if __name__ == "__main__":
    wf_id = build_workflow()
    print(wf_id)
