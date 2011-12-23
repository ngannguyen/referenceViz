#!/usr/bin/env python

"""
Input: pileup file
Output: coverage distribution
"""

import os, sys, re

cov2count = {}
total = 0
for line in sys.stdin:
    items = line.strip().split()
    cov = int(items[3])
    total += 1 
    if cov not in cov2count:
        cov2count[cov] = 1
    else:
        cov2count[cov] += 1

for cov in cov2count:
    freq = cov2count[cov]*100.0/total
    sys.stdout.write("%d\t%f\n" %(cov, freq))
