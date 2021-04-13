#!/usr/bin/env python

import os
import sys


sample_files_list = sys.argv[1]
count_files_lists = sys.argv[2]
out_counts_file = sys.argv[3]
out_samples_file = sys.argv[4]
OUTCOUNTS = open(out_counts_file, 'w')
OUTSAMPLES = open(out_samples_file, 'w')
ERRORS = open("errors.txt", 'w')

SAMPLEFILES = open(sample_files_list)
COUNTFILES = open(count_files_lists)

sample_header = ''
sample_order = 0
samples = {}
sample_lines = {}
feature_list = {}

print("Looping SAMPLEFILES")
for sample_file in SAMPLEFILES:

    sample_file = sample_file.rstrip("\n")
    SAMPLEFILE = open(sample_file)

    sample_header = ''
    line_counter = 0
    for line in SAMPLEFILE:
        line = line.rstrip("\n")
        if line[0] == "#":
            continue
        if line_counter == 0:
            tmp_sample_header = line
            if len(sample_header) > 0:
                if tmp_sample_header != sample_header:
                    ERRORS.write("SAMPLE HEADERS DO NOT MATCH: " + sample_header + " vs " + tmp_sample_header)
                    sys.exit(0)
            else:
                sample_header = tmp_sample_header
            line_counter += 1
            continue
        line_counter += 1
        line_data = line.split("\t")
        sample_name = line_data[0]
        if sample_name not in samples:
            samples[sample_name] = sample_order
            sample_lines[sample_name] = line
        sample_order += 1

counter_samples = {}
got_feature_list = False
print("done looping sample files\n")
print("Samples:")
print(samples)
print("")
count_results = {}

for count_file in COUNTFILES:
    count_file = count_file.rstrip("\n")
    COUNTFILE = open(count_file)
    count_header = COUNTFILE.readline()
    count_header = count_header.rstrip("\n")
    count_header_list = count_header.split("\t")
    print(count_header)
    for feature_num, count_line in enumerate(COUNTFILE):
        count_line = count_line.rstrip("\n")
        count_data = count_line.split("\t")
        feature = count_data[0]
        if feature[0:2] == "__":
            continue
        # if we arleady got a feature list .. check if same order
        if got_feature_list:
            old_feature = feature_list[feature_num]
            # features should be the same
            if feature != old_feature:
                print("FEATURE LISTS BETWEEN COUNT FILES NOT THE SAME")
                ERRORS.write("FEATURE LISTS BETWEEN COUNT FILES NOT THE SAME")
                sys.exit()
        else:
            # if no old feature list put feature in new list
            feature_list[feature_num] = feature
            count_results[feature] = {}

        for sc in range(1, len(count_data)):
            sample = count_header_list[sc]
            if sample not in samples:
                print("SAMPLE " + sample + " not in a sample sheet\nDoes the file have a header?")
                ERRORS.write("SAMPLE " + sample + " not in a sample sheet\nDoes the file have a header?")
                sys.exit()

            count = count_data[sc]
            if sample in count_results[feature]:
                if count != count_results[feature][sample]:
                    print("TWO SAMPLES WITH SAME NAME HAVE DIFFERENT COUNT VALUES")
                    ERRORS.write("TWO SAMPLES WITH SAME NAME HAVE DIFFERENT COUNT VALUES")
                    sys.exit()
            else:
                count_results[feature][sample] = count
    got_feature_list = True
print("done getting features")
samples_ordered = []
for key, value in sorted(samples.iteritems(), key=lambda (k, v): (v, k)):
    samples_ordered.append(key)
header_line = count_header_list[0]
for s in samples_ordered:
    header_line += "\t" + s

print("SAMPLES", samples_ordered)
OUTCOUNTS.write(header_line + "\n")
for feature_c in sorted(feature_list.keys()):
    feature = feature_list[feature_c]
    out_line = feature
    for sample in samples_ordered:
        if sample not in count_results[feature]:
            ERRORS.write("COUNTS NOT AVAIABLE FOR SAMPLE " + sample)
            sys.exit()
        count = count_results[feature][sample]
        out_line += "\t" + count
    OUTCOUNTS.write(out_line + "\n")

OUTSAMPLES.write(sample_header + "\n")
for s in samples_ordered:
    sample_line = sample_lines[s]
    OUTSAMPLES.write(sample_line + "\n")
