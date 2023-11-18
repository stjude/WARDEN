#!/usr/bin/env python3

import sys
import os
import re
from decimal import *

sample_list = sys.argv[1]
SAMPLELIST = open(sample_list)
header = SAMPLELIST.readline()

out_file = "alignment_statistics.txt"
OUT = open(out_file, "w")

cols = (
    "Sample",
    "totalReads",
    "totalQCPassed",
    "totalQCNotPassed",
    "Duplicate",
    "DuplicatePerc",
    "MappedReads",
    "PercentMapped",
    "Paired",
)

OUT.write("\t".join(cols) + "\n")


for line in SAMPLELIST:
    line = line.rstrip("\n")
    line_data = line.split()
    sample_name = line_data[0]
    flagstat_file = "FLAGSTAT/%s.Aligned.flagstat.txt" % sample_name
    FLAG = open(flagstat_file)
    flag_lines = FLAG.read()
    m1 = re.search("(\d+)\s+\+\s+(\d+)\s+in total", flag_lines)
    totalQCPassed = int(m1.group(1))
    totalQCNotPassed = int(m1.group(2))
    totalReads = totalQCPassed + totalQCNotPassed
    m2 = re.search("(\d+)\s+\+\s+(\d+)\s+secondary", flag_lines)
    m3 = re.search("(\d+)\s+\+\s+(\d+)\s+duplicates", flag_lines)
    duplicates = int(m3.group(1))

    m4 = re.search("(\d+)\s+\+\s+(\d+)\s+mapped \((\S+)", flag_lines)
    mapped = m4.group(1)
    duplicate_perc = Decimal(duplicates) / Decimal(mapped)
    duplicate_perc = duplicate_perc * 100
    per_mapped = m4.group(3)
    m5 = re.search("(\d+)\s+\+\s+(\d+)\s+paired in sequencing", flag_lines)
    paired = m5.group(1)
    out_data = (
        sample_name,
        str(totalReads),
        str(totalQCPassed),
        str(totalQCNotPassed),
        str(duplicates),
        str(duplicate_perc),
        mapped,
        per_mapped,
        paired,
    )
    OUT.write("\t".join(out_data) + "\n")
