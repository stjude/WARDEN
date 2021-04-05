#!/usr/bin/env python
# warden_python 1.0.0


import sys
import dxpy
import re
import json

parameter_file = sys.argv[1]
fastq_samples_file = sys.argv[2]
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
FQFILES = open(fastq_samples_file)
for line in FQFILES:
    line = line.rstrip("\n")
    sample_name, for_fastq, rev_fastq = line.split("\t")
    samples[sample_name] = [for_fastq, rev_fastq]


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

    fastqc_applet = dxpy.search.find_one_data_object(classname="applet",
                                                     name=app_names["fastqc"],
                                                     state="closed",
                                                     return_handler=True)
    star_applet = dxpy.search.find_one_data_object(classname="applet",
                                                   name=app_names["star"],
                                                   state="closed",
                                                   return_handler=True)
    combine_sj_tab_applet = dxpy.search.find_one_data_object(classname="applet",
                                                             name=app_names["combine_sj_out"],
                                                             state="closed",
                                                             return_handler=True)
    sort_bam_applet = dxpy.search.find_one_data_object(classname="applet",
                                                       name=app_names["sort_bam"],
                                                       state="closed",
                                                       return_handler=True)
    htseq_applet = dxpy.search.find_one_data_object(classname="applet",
                                                    name=app_names["htseq"],
                                                    state="closed",
                                                    return_handler=True)
    genome_cov_applet = dxpy.search.find_one_data_object(classname="applet",
                                                         name=app_names["genome_coverage"],
                                                         state="closed",
                                                         return_handler=True)
    bigwig_applet = dxpy.search.find_one_data_object(classname="applet",
                                                     name=app_names["bigwig"],
                                                     state="closed",
                                                     return_handler=True)
    combine_counts_applet = dxpy.search.find_one_data_object(classname="applet",
                                                             name=app_names["combine_counts"],
                                                             state="closed",
                                                             return_handler=True)
    combine_flagstat_applet = dxpy.search.find_one_data_object(classname="applet",
                                                               name=app_names["combine_flagstat"],
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
    bw_viewer_applet = dxpy.search.find_one_data_object(classname="applet",
                                                        name=app_names["bw_viewer"],
                                                        state="closed",
                                                        return_handler=True)

    sample_num = 0
    htseq_results = []
    bigwig_files = []
    flagstat_files_arr = []
    index_project, index_id = parameters["index_file"].split(":")
    gtf_project, gtf_id = parameters["gtf_file"].split(":")
    genome_length_project, genome_length_id = parameters["genome_sizes_file"].split(":")
    gene_length_project, gene_length_id = parameters["gene_length_file"].split(":")
    fpkm_results = []
    fpkm_log2_results = []
    sj_out_files_arr = []

    star_alignment_opts = {
        "outSAMunmapped": parameters["outSAMunmapped"],
        "outSAMattributes": parameters["outSAMattributes"],
        "outFilterMultimapNmax": int(parameters["outFilterMultimapNmax"]),
        "outFilterMismatchNmax": int(parameters["outFilterMismatchNmax"]),
        "alignIntronMax": int(parameters["alignIntronMax"]),
        "outSAMstrandField": parameters["outSAMstrandField"],
        "chimSegmentMin": int(parameters["chimSegmentMin"]),
        "sjdbOverhang": int(parameters["sjdbOverhang"]),
        "chimJunctionOverhangMin": int(parameters["chimJunctionOverhangMin"]),
        "subsample_target": int(parameters["STAR_subsample_n_reads"]),
    }

    if parameters["two_pass_alignment"] == 'true':
        for sample_name in samples:
            forward_id = samples[sample_name][0]
            forward_link = dxpy.dxlink(forward_id)
            align_input = {}
            align_input["first_pass"] = True
            align_input.update(star_alignment_opts)
            align_input["read_file1"] = forward_link
            star_instance = parameters["star_instance"]
            align_input["mark_duplicates"] = False
            align_input["generate_transcriptome_BAM"] = False
            align_input["star_index_archive"] = dxpy.dxlink({"project": index_project, "id": index_id})
            if parameters["sjdbFileChrStartEnd"] != "null" and parameters["sjdbFileChrStartEnd"] != '':
                sjdbFileChrStartEnd_project, sjdbFileChrStartEnd_id = parameters["sjdbFileChrStartEnd"].split(":")
                align_input["sjdbFileChrStartEnd"] = dxpy.dxlink({"project": sjdbFileChrStartEnd_project, "id": sjdbFileChrStartEnd_id})
            if parameters["indexed_with_gtf"] != "true":
                align_input["transcriptome_gtf"] = dxpy.dxlink({"project": gtf_project, "id": gtf_id})
            align_input["output_prefix"] = sample_name
            align_stage_id = ""
            if samples[sample_name][1] != "-":
                reverse_id = samples[sample_name][1]
                reverse_link = dxpy.dxlink(reverse_id)
                align_input["read_file2"] = reverse_link
                first_align_stage_id = wf.add_stage(star_applet, stage_input=align_input, instance_type=star_instance, folder="ALIGN_Pass1", name=sample_name + ":ALIGN_Pass1")
            else:
                first_align_stage_id = wf.add_stage(star_applet, stage_input=align_input, instance_type=star_instance, folder="ALIGN_Pass1", name=sample_name + ":ALIGN_Pass1")

            sj_out_files_arr.append(dxpy.dxlink({"stage": first_align_stage_id, "outputField": "sj_tab_out"}))
        combine_sj_out_input = {"sj_out_files": sj_out_files_arr}
        combine_sj_out_stage_pass1_id = wf.add_stage(combine_sj_tab_applet, stage_input=combine_sj_out_input, instance_type="azure:mem2_ssd1_x1", name="COMBINE SJ OUT PASS1", folder="COMBINED_JUNCTIONS_PASS1")
        parameters["pass1_sj_out"] = dxpy.dxlink({"stage": combine_sj_out_stage_pass1_id, "outputField": "combined_sj_out"})

    for sample_name in samples:
        forward_id = samples[sample_name][0]
        forward_link = dxpy.dxlink(forward_id)
        if parameters["run_FastQC"] == 'true':
            forward_input = {"fastq_input": forward_link}
            fq_stage_id = wf.add_stage(fastqc_applet, stage_input=forward_input, instance_type="azure:mem2_ssd1_x2", folder="FASTQC", name=sample_name + ":Forward FASTQC")
            wf_outputs += [
                {
                    "name": sample_name + "_forward_fastqc_html",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": fq_stage_id,
                            "outputField": "html_file"
                        }
                    }
                },
                {
                    "name": sample_name + "_forward_fastqc_zip",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": fq_stage_id,
                            "outputField": "zip_file"
                        }
                    }
                }
            ]

        align_input = {}
        align_input.update(star_alignment_opts)
        mark_duplicates = parameters["mark_duplicates"]
        if mark_duplicates == "false":
            mark_duplicates = False
        else:
            mark_duplicates = True
        align_input["mark_duplicates"] = mark_duplicates
        if parameters["generate_transcriptome_BAM"] == "true":
            align_input["generate_transcriptome_BAM"] = True
        else:
            align_input["generate_transcriptome_BAM"] = False
        star_instance = parameters["star_instance"]
        align_input["read_file1"] = forward_link
        align_input["star_index_archive"] = dxpy.dxlink({"project": index_project, "id": index_id})
        if parameters["indexed_with_gtf"] != "true":
            align_input["transcriptome_gtf"] = dxpy.dxlink({"project": gtf_project, "id": gtf_id})
        if "pass1_sj_out" in parameters:
            align_input["sjdbFileChrStartEnd"] = parameters["pass1_sj_out"]
        elif parameters["sjdbFileChrStartEnd"] != "null" and parameters["sjdbFileChrStartEnd"] != '':
            sjdbFileChrStartEnd_project, sjdbFileChrStartEnd_id = parameters["sjdbFileChrStartEnd"].split(":")
            align_input["sjdbFileChrStartEnd"] = dxpy.dxlink({"project": sjdbFileChrStartEnd_project, "id": sjdbFileChrStartEnd_id})
        align_input["output_prefix"] = sample_name
        align_stage_id = ""
        if samples[sample_name][1] != "-":
            reverse_id = samples[sample_name][1]
            reverse_link = dxpy.dxlink(reverse_id)
            if parameters["run_FastQC"] == 'true':
                rev_input = {"fastq_input": reverse_link}
                rev_fq_stage_id = wf.add_stage(fastqc_applet, stage_input=rev_input, instance_type="azure:mem2_ssd1_x2", folder="FASTQC", name=sample_name + ":Reverse FASTQC")
                wf_outputs += [
                    {
                        "name": sample_name + "_reverse_fastqc_html",
                        "class": "file",
                        "outputSource": {
                            "$dnanexus_link": {
                                "stage": rev_fq_stage_id,
                                "outputField": "html_file"
                            }
                        }
                    },
                    {
                        "name": sample_name + "_reverse_fastqc_zip",
                        "class": "file",
                        "outputSource": {
                            "$dnanexus_link": {
                                "stage": rev_fq_stage_id,
                                "outputField": "zip_file"
                            }
                        }
                    }
                ]

            align_input["read_file2"] = reverse_link
            align_stage_id = wf.add_stage(star_applet, stage_input=align_input, instance_type=star_instance, folder="STAR", name=sample_name + ":ALIGN")
        else:
            align_stage_id = wf.add_stage(star_applet, stage_input=align_input, instance_type=star_instance, folder="STAR", name=sample_name + ":ALIGN")
        flagstat_files_arr.append(dxpy.dxlink({"stage": align_stage_id, "outputField": "flagstat_out"}))
        sj_out_files_arr.append(dxpy.dxlink({"stage": align_stage_id, "outputField": "sj_tab_out"}))
        wf_outputs += [
            {
                "name": sample_name + "_star_bam",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "sorted_by_coord_bam"
                    }
                }
            },
            {
                "name": sample_name + "_star_log",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "log_final_out"
                    }
                }
            },
            {
                "name": sample_name + "_flagstat",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "flagstat_out"
                    }
                }
            },
            {
                "name": sample_name + "_star_splice_junctions",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "sj_tab_out"
                    }
                }
            },
            {
                "name": sample_name + "_star_chimeric_bam",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "chimeric_bam"
                    }
                }
            },
            {
                "name": sample_name + "_star_chimeric_junction",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": align_stage_id,
                        "outputField": "chimeric_junction"
                    }
                }
            },
        ]
        if parameters["generate_transcriptome_BAM"] == "true":
            wf_outputs += [
                {
                    "name": sample_name + "_star_transcriptome_bam",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": align_stage_id,
                            "outputField": "to_transcriptome_bam"
                        }
                    }
                },
            ]

        if parameters["generate_name_sorted_BAM"] == "true":
            sort_input = {"input_bam": dxpy.dxlink({"stage": align_stage_id, "outputField": "sorted_by_coord_bam"})}
            sort_stage_id = wf.add_stage(sort_bam_applet, stage_input=sort_input, instance_type="azure:mem2_ssd1_x2", name=sample_name + ":NAME SORT BAM", folder="STAR")
            wf_outputs += [
                {
                    "name": sample_name + "_name_sorted_bam",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": sort_stage_id,
                            "outputField": "output_bam"
                        }
                    }
                },
            ]
            htseq_input = {"input_bam": dxpy.dxlink({"stage": sort_stage_id, "outputField": "output_bam"})}
            htseq_input["order"] = "name"
        else:
            htseq_input = {"input_bam": dxpy.dxlink({"stage": align_stage_id, "outputField": "sorted_by_coord_bam"})}
            htseq_input["order"] = "pos"

        htseq_input["annotation_file"] = dxpy.dxlink({"project": gtf_project, "id": gtf_id})
        htseq_input["gene_length_file"] = dxpy.dxlink({"project": gene_length_project, "id": gene_length_id})
        htseq_input["prefix"] = sample_name
        htseq_input["strand"] = parameters["strandedness"]
        htseq_input["feature_type"] = parameters["feature_type"]
        htseq_input["id_attribute"] = parameters["id_attribute"]
        htseq_input["mode"] = parameters["mode"]
        htseq_input["nonunique"] = parameters["nonunique"]
        htseq_input["secondary_alignments"] = parameters["secondary_alignments"]
        htseq_input["supplementary_alignments"] = parameters["supplementary_alignments"]
        htseq_stage_id = wf.add_stage(htseq_applet, stage_input=htseq_input, instance_type=parameters["htseq_instance"], name=sample_name + ":HTSEQ COUNT", folder="HTSEQ")
        htseq_results.append(dxpy.dxlink({"stage": htseq_stage_id, "outputField": "htseq_counts"}))
        wf_outputs += [
            {
                "name": sample_name + "_htseqcounts",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": htseq_stage_id,
                        "outputField": "htseq_counts"
                    }
                }
            },
        ]
        if parameters["id_attribute"] == "gene_name":
            fpkm_results.append((dxpy.dxlink({"stage": htseq_stage_id, "outputField": "fpkm"})))
            fpkm_log2_results.append((dxpy.dxlink({"stage": htseq_stage_id, "outputField": "fpkm_log2"})))
            wf_outputs += [
                {
                    "name": sample_name + "_fpkm",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": htseq_stage_id,
                            "outputField": "fpkm"
                        }
                    }
                },
                {
                    "name": sample_name + "_fpkm_log2",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": htseq_stage_id,
                            "outputField": "fpkm_log2"
                        }
                    }
                },
            ]

        if parameters["run_coverage"] == 'true':
            gcb_input = {}
            gcb_input["input_bam"] = dxpy.dxlink({"stage": align_stage_id, "outputField": "sorted_by_coord_bam"})
            gcb_input["genome_sizes_file"] = dxpy.dxlink({"project": genome_length_project, "id": genome_length_id})
            gcb_input["strandedness"] = parameters["strandedness"]
            gcb_input["output_prefix"] = sample_name
            gcb_stage_id = wf.add_stage(genome_cov_applet, stage_input=gcb_input, instance_type="azure:mem3_ssd1_x8", name=sample_name + ":COVERAGE", folder="COVERAGE")

            bg2bw_all_input = {}
            bg2bw_all_input["bedgraph_file"] = dxpy.dxlink({"stage": gcb_stage_id, "outputField": "all_coverage_file"})
            bg2bw_all_input["genome_sizes_file"] = dxpy.dxlink({"project": genome_length_project, "id": genome_length_id})
            bg2bw_all_input["output_prefix"] = sample_name
            bg2bw_all_stage_id = wf.add_stage(bigwig_applet, stage_input=bg2bw_all_input, instance_type="azure:mem2_ssd1_x4", name=sample_name + ":BED To BW-ALL", folder="BIGWIG")
            bigwig_files.append(dxpy.dxlink({"stage": bg2bw_all_stage_id, "outputField": "bigwig"}))
            wf_outputs += [
                {
                    "name": sample_name + "_all_bigwig",
                    "class": "file",
                    "outputSource": {
                        "$dnanexus_link": {
                            "stage": bg2bw_all_stage_id,
                            "outputField": "bigwig"
                        }
                    }
                },
            ]

            if parameters["strandedness"] != "no":
                bg2bw_pos_input = {}
                bg2bw_pos_input["bedgraph_file"] = dxpy.dxlink({"stage": gcb_stage_id, "outputField": "pos_coverage_file"})
                bg2bw_pos_input["genome_sizes_file"] = dxpy.dxlink({"project": genome_length_project, "id": genome_length_id})
                bg2bw_pos_input["output_prefix"] = sample_name
                bg2bw_pos_stage_id = wf.add_stage(bigwig_applet, stage_input=bg2bw_pos_input, instance_type="azure:mem2_ssd1_x4", name=sample_name + ":BED To BW-POS", folder="BIGWIG")
                wf_outputs += [
                    {
                        "name": sample_name + "_pos_bigwig",
                        "class": "file",
                        "outputSource": {
                            "$dnanexus_link": {
                                "stage": bg2bw_pos_stage_id,
                                "outputField": "bigwig"
                            }
                        }
                    },
                ]

                bg2bw_neg_input = {}
                bg2bw_neg_input["bedgraph_file"] = dxpy.dxlink({"stage": gcb_stage_id, "outputField": "neg_coverage_file"})
                bg2bw_neg_input["genome_sizes_file"] = dxpy.dxlink({"project": genome_length_project, "id": genome_length_id})
                bg2bw_neg_input["output_prefix"] = sample_name
                bg2bw_neg_stage_id = wf.add_stage(bigwig_applet, stage_input=bg2bw_neg_input, instance_type="azure:mem2_ssd1_x4", name=sample_name + ":BED To BW-NEG", folder="BIGWIG")
                wf_outputs += [
                    {
                        "name": sample_name + "_neg_bigwig",
                        "class": "file",
                        "outputSource": {
                            "$dnanexus_link": {
                                "stage": bg2bw_neg_stage_id,
                                "outputField": "bigwig"
                            }
                        }
                    },
                ]

                bigwig_files.append(dxpy.dxlink({"stage": bg2bw_pos_stage_id, "outputField": "bigwig"}))
                bigwig_files.append(dxpy.dxlink({"stage": bg2bw_neg_stage_id, "outputField": "bigwig"}))
        sample_num += 1

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
    if parameters["id_attribute"] == "gene_name":
        combine_fpkm_input = {"count_files": fpkm_results, "name_value": "fpkm", "sample_files": [dxpy.dxlink(final_sample_list_id)]}
        combine_fpkm_stage_id = wf.add_stage(combine_counts_applet, stage_input=combine_fpkm_input, instance_type="azure:mem2_ssd1_x1", name="COMBINE FPKM")
        combine_fpkm_log2_input = {"count_files": fpkm_log2_results, "name_value": "fpkm.log2", "sample_files": [dxpy.dxlink(final_sample_list_id)]}
        combine_fpkm_log2_stage_id = wf.add_stage(combine_counts_applet, stage_input=combine_fpkm_log2_input, instance_type="azure:mem2_ssd1_x1", name="COMBINE FPKMlog2")
        wf_outputs += [
            {
                "name": "combined_fpkm",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": combine_fpkm_stage_id,
                        "outputField": "count_file"
                    }
                }
            },
            {
                "name": "combined_fpkm_log2",
                "class": "file",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": combine_fpkm_log2_stage_id,
                        "outputField": "count_file"
                    }
                }
            },
        ]

    combine_flagstat_input = {"flagstat_files": flagstat_files_arr, "sample_list": dxpy.dxlink(final_sample_list_id)}
    combine_flagstat_stage_id = wf.add_stage(combine_flagstat_applet, stage_input=combine_flagstat_input, instance_type="azure:mem2_ssd1_x1", name="COMBINE FLAGSTAT", folder="STAR")
    wf_outputs += [
        {
            "name": "combined_flagstat",
            "class": "file",
            "outputSource": {
                "$dnanexus_link": {
                    "stage": combine_flagstat_stage_id,
                    "outputField": "combined_flagstat"
                }
            }
        },
    ]

    if parameters["BW_VIEWER"] != "None" and parameters["run_coverage"] == 'true':
        bw_project, bw_file = parameters["BW_VIEWER"].split(":")
        viewer_link = dxpy.dxlink({"project": bw_project, "id": bw_file})
        bw_viewer_input = {"viewer": viewer_link, "bigwig_files": bigwig_files}
        bw_viewer_stage_id = wf.add_stage(bw_viewer_applet, stage_input=bw_viewer_input, instance_type="azure:mem2_ssd1_x1", name="BIGWIG_VIEWER", folder="BIGWIG")
        wf_outputs += [
            {
                "name": "bw_viewer",
                "class": "record",
                "outputSource": {
                    "$dnanexus_link": {
                        "stage": bw_viewer_stage_id,
                        "outputField": "viewer_bookmark"
                    }
                }
            },
        ]

    if parameters["limma_DE_viewer"] != "None":
        limma_viewer_project, limma_viewer_file = parameters["limma_DE_viewer"].split(":")
        limma_viewer_link = dxpy.dxlink({"project": limma_viewer_project, "id": limma_viewer_file})

    if parameters["run_limma"] == 'true' and parameters["limma_runnable"] == "true":
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
    if parameters["run_simple_dif_ex"] == 'true':
        simple_DE_input = {"input_count_file": dxpy.dxlink({"stage": combine_counts_stage_id, "outputField": "count_file"}), "sample_list_file": dxpy.dxlink(final_sample_list_id), "contrasts_file": dxpy.dxlink(comparisons_all_id), "difex_viewer": limma_viewer_link}
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
