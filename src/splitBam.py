#!/usr/bin/env python

import os, sys, re
import pysam

if len(sys.argv) < 3:
    sys.stderr.write("Usage: splitBam.py <numLinePerFile> <prefix> < input.bam\n")
    sys.exit(1)

linePerFile = int(sys.argv[1])
prefix = sys.argv[2]

f = pysam.Samfile("-", "rb")

suffix = 0
filename = "%s-%d.bam" %(prefix, suffix)
outfile = pysam.Samfile( filename, "wb", template=f )

count = 0
for read in f:
    if count < linePerFile:
        count += 1
        outfile.write(read)
    else:
        outfile.close()
        suffix += 1
        filename = "%s-%d.bam" %(prefix, suffix)
        outfile = pysam.Samfile( filename, "wb", template = f )
        outfile.write(read)
        count = 1

outfile.close()
f.close()
