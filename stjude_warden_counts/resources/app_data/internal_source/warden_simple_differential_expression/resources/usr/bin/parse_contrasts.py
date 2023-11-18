#!/usr/bin/env python3
import sys


contrasts = sys.argv[1]
contrasts = "".join(contrasts.split())
comp_array = contrasts.split(",")
out_comps = []
for c in comp_array:
    c1, c2 = c.split("-")
    out_comps.append(c1 + "_vs_" + c2 + "=" + c1 + "-" + c2)
print(",".join(out_comps))
