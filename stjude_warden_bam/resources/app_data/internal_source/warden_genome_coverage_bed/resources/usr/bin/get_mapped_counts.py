#!/usr/bin/env python3
import sys
import re

flagstat = sys.argv[1]
IN = open(flagstat)
for line in IN:
    m = re.search("(\d+)\s+\+\s+\d+\s+mapped", line)
    if m:
        print(m.group(1))
        sys.exit()
print("0")
