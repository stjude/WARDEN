#!/usr/bin/env python

import sys
from decimal import *
import math

count_file = sys.argv[1]
gene_length_file = sys.argv[2]
out_file = sys.argv[3]
out_file_log2 = sys.argv[4]
OUT = open(out_file, "w")
OUTLOG2 = open(out_file_log2, "w")
COUNTS = open(count_file)
total_hits = 0
count_header = COUNTS.readline()
for line in COUNTS:
    line = line.rstrip("\n")
    line = line.rstrip("\r")
    g, count = line.split("\t")
    if g == "__not_aligned":
        continue
    total_hits += int(count)
COUNTS.close()
print("total_hits", total_hits)
GL = open(gene_length_file)
gene_lengths = {}
for line in GL:
    line = line.rstrip("\n")
    line = line.rstrip("\r")
    gene, gLength = line.split("\t")
    gene_lengths[gene] = gLength

scaling_factor = Decimal(total_hits) / Decimal(1000000)
print("Scaling", scaling_factor)
COUNTS = open(count_file)
total_hits = 0
count_header = COUNTS.readline()
OUT.write(count_header)
OUTLOG2.write(count_header)
for line in COUNTS:
    line = line.rstrip("\n")
    line = line.rstrip("\r")
    gene, count = line.split("\t")
    if gene[0] == "_":
        continue
    rpm = Decimal(count) / scaling_factor
    rpmlog2 = (Decimal(count) + Decimal(0.5)) / scaling_factor
    gene_length = gene_lengths[gene]
    if gene_length == 0:
        continue
    fpkm = rpm / (Decimal(gene_length) / Decimal(1000))
    fpkmlog2 = rpmlog2 / (Decimal(gene_length) / Decimal(1000))
    log2fpkm = math.log(fpkmlog2, 2)
    OUT.write("%s\t%.3f\n" % (gene, fpkm))
    OUTLOG2.write("%s\t%.3f\n" % (gene, log2fpkm))
