#!/usr/bin/env python3

import sys
import os
import re
import subprocess

standalone = sys.argv[1]
sample_list = sys.argv[2]  # samplelist file
fastq_list = sys.argv[3]  # path to fastqs
fastq_object_list = sys.argv[4]  # object id

SAMPLELIST = open(sample_list)
FASTQLIST = open(fastq_list)
FASTQOBJECTLIST = open(fastq_object_list)


def safe_name(name):
    if not name[0].isalpha():
        name = "x" + name
    outname = ""
    for char in name:
        if not char.isalnum() and char != "_":
            outname += "_"
        else:
            outname += char
    return outname


def process_comparisons(comparison_line, comparison_list):
    comparison_split = comparison_line.strip().split("=")
    if len(comparison_split) != 2:
        PROCESSERRORS.write(
            "Comparison line should be of following format: `#COMPARISONS=PHENOTYPE1-PHENOTYPE2` or `#COMPARISONS=PHENOTYPE1-PHENOTYPE2,PHENOTYPE1-PHENOTYPE3`"
        )
        PROCESSERRORS.close()
        sys.exit()
    _, comparison_list_string = comparison_split
    comparisons = comparison_list_string.split(",")
    for comparison in comparisons:
        comparison_split = comparison.split("-")
        if len(comparison_split) != 2:
            PROCESSERRORS.write(
                "Comparison line should be of following format: `#COMPARISONS=PHENOTYPE1-PHENOTYPE2` or `#COMPARISONS=PHENOTYPE1-PHENOTYPE2,PHENOTYPE1-PHENOTYPE3`"
            )
            PROCESSERRORS.close()
            sys.exit()
        comparison_list.append(comparison_split)


fastq_objects = {}
for line in FASTQLIST:
    fastq_path = line.rstrip("\n")
    fastq_file = fastq_path.split("/")[-1]
    fastq_file_object = FASTQOBJECTLIST.readline().rstrip("\n")
    fastq_objects[fastq_file] = fastq_file_object

FASTQSAMPLELIST = open("bam_file_list.txt", "w")
FINALLIST = open("cleaned_sample_list.txt", "w")
COMPARISONLIMMA = open("comparisons_limma.txt", "w")
COMPARISON_ALL = open("comparisons_simple_DE.txt", "w")
PROCESSERRORS = open("process_files_errors.txt", "w")

comparison_list = []
fastqs = []
phenotypes = {}
line_counter = 0
for line in SAMPLELIST:
    line = line.rstrip("\n")
    if len("".join(line.split())) == 0:
        continue
    if line[0] == "#":
        if (
            line[0:13].upper() == "#COMPARISONS="
            or line[0:12].upper() == "#COMPARISON="
        ):
            process_comparisons(line, comparison_list)
        continue

    if line_counter == 0:
        header = line
        header_split = header.split()
        header_count = len(header_split)
        if header_count < 3:
            PROCESSERRORS.write(
                "The sample list file does not have at least three columns:"
                + header
                + "\n"
            )
            PROCESSERRORS.close()
            sys.exit()
        FINALLIST.write("\t".join(header_split[:3]) + "\n")
        line_counter += 1
        continue

    line_counter += 1
    line_data = line.split()
    num_cols = len(line_data)
    if num_cols == 0:
        continue
    if num_cols < 3:
        PROCESSERRORS.write(
            "Line "
            + line_counter
            + " of the sample list does not have at least three columns\n"
        )
        PROCESSERRORS.close()
        sys.exit()

    sample_name = safe_name(line_data[0])
    read_file1 = line_data[2]

    if read_file1 in fastqs:
        PROCESSERRORS.write(read_file1 + " listed twice\n")
        PROCESSERRORS.close()
        sys.exit()

    fastqs.append(read_file1)
    phenotype = safe_name(line_data[1])
    if phenotype not in phenotypes:
        phenotypes[phenotype] = 0
    phenotypes[phenotype] += 1

    if read_file1 not in fastq_objects:
        PROCESSERRORS.write(
            read_file1 + " does not exist but is listed in the sample list file\n"
        )
        PROCESSERRORS.close()
        sys.exit()

    object1 = fastq_objects[read_file1]
    FINALLIST.write("\t".join([sample_name, phenotype, read_file1]) + "\n")
    FASTQSAMPLELIST.write(sample_name + "\t" + object1 + "\n")

if standalone == "true":
    limma_comparisons = []
    all_comparisons = []
    for c in comparison_list:
        c1, c2 = c
        c1 = safe_name(c1)
        c2 = safe_name(c2)
        if c1 not in phenotypes:
            PROCESSERRORS.write(c1 + " from comparisons not found in a phenotype\n")
            PROCESSERRORS.close()
            sys.exit()
        if c2 not in phenotypes:
            PROCESSERRORS.write(c2 + " from comparisons not found in a phenotype\n")
            PROCESSERRORS.close()
            sys.exit()
        if c1 == c2:
            PROCESSERRORS.write(
                "Can't compare "
                + c1
                + " and "
                + c2
                + " because they're the same phenotype\n"
            )
            PROCESSERRORS.close()
            sys.exit()
        c1Count = phenotypes[c1]
        c2Count = phenotypes[c2]
        if c1Count >= 2 and c2Count >= 2:
            limma_comparisons.append(c1 + "-" + c2)
        all_comparisons.append(c1 + "-" + c2)

    COMPARISONLIMMA.write(",".join(limma_comparisons) + "\n")
    COMPARISON_ALL.write(",".join(all_comparisons) + "\n")

for filename in fastq_objects.keys():
    if filename not in fastqs:
        PROCESSERRORS.write("BAM file not found in sample file:" + filename + "\n")
        PROCESSERRORS.close()
        sys.exit()
