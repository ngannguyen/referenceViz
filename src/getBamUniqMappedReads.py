#!/usr/bin/env python

import os, sys, re
import pysam

f = pysam.Samfile("-", "rb") #Input bam from stdin
outfile = pysam.Samfile("-", "wb", template=f) #Output bam to stdin

#count = 0
for r in f:
    if r.tags != None and 'XA' in [ tag for (tag, vals) in r.tags ]:
        continue
    outfile.write( r )
    #count += 1

outfile.close()
f.close()

