#!/usr/bin/env python

import os, sys, re
import pysam

if len(sys.argv) < 2:
    sys.stderr.write("Usage: splitBam.py <prefix> < input.bam\n")
    sys.exit(1)

prefix = sys.argv[1]

f = pysam.Samfile("-", "rb")

pairedfilename = "%s-paired.bam" %(prefix)
pairedoutfile = pysam.Samfile( pairedfilename, "wb", template=f )
singlefilename = "%s-single.bam" %(prefix)
singleoutfile = pysam.Samfile( singlefilename, "wb", template=f )

currread = None
for read in f:
    if currread != None and read.qname == currread.qname:
        if (currread.is_read1 and read.is_read2) or (currread.is_read2 and read.is_read1):
            pairedoutfile.write(read)
            pairedoutfile.write(currread)
            currread = None
        else:
            sys.stderr.write("Not mate reads but have same name!\n")
            sys.exit(1)
    elif currread != None:
        singleoutfile.write(currread)
        currread = read
    else:
        currread = read

if currread != None:
    singleoutfile.write(currread)

pairedoutfile.close()
singleoutfile.close()
f.close()
