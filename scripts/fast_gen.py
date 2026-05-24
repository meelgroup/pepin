#!/usr/bin/env python3
# Quick DNF generator (faster than the scaffolded one for large M).
# Usage: fast_gen.py <nvars> <nclauses> <clsize> <out.dnf> [seed]
import sys, random
nvars = int(sys.argv[1])
ncls  = int(sys.argv[2])
ksize = int(sys.argv[3])
out   = sys.argv[4]
seed  = int(sys.argv[5]) if len(sys.argv) > 5 else 1
random.seed(seed)
with open(out, 'w') as f:
    f.write(f"p dnf {nvars} {ncls}\n")
    vars_list = list(range(1, nvars+1))
    for _ in range(ncls):
        sv = random.sample(vars_list, ksize)
        f.write(" ".join(str(v) for v in sv) + " 0\n")
