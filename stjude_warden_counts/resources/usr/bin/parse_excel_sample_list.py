#!/usr/bin/env python

import openpyxl
import csv
import sys

wb = openpyxl.load_workbook(sys.argv[1])
sh = wb.active
c = csv.writer(
    sys.stdout, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="", quotechar=""
)
for r in sh.rows:
    c.writerow([cell.value for cell in r])
