#!/usr/bin/env python

import os, sys, re, copy
import pysam

f = pysam.Samfile("-", "rb") #Input bam from stdin
chr = 6
start = 28477754

header = { 'HD': {'VN': '1.0'},
           'SQ': [{'LN': 171115067, 'SN': '6'}]}
out = pysam.Samfile("-", "wb", header=header) #Output bam to stdout

for r in f:
    r.rname = 0
    r.pos += start
    r.mpos += start
    out.write(r)
    #hg19.chr6.171115067.28477754.4970600.1
    #sys.stderr.write("%s\n" % r.rname )
    #items = name.split('.')
    #newread = copy.copy(r)
    #a = pysam.AlignedRead()
    #a.qname = r.qname
    #a.seq= r.seq
    #a.flag = r.flag
    #a.rname = chr
    #a.pos = r.pos + start
    #a.mapq = r.mapq
    #a.cigar = r.cigar
    #a.mrnm = r.mrnm
    #a.mpos= r.mpos + start
    #a.isize= r.isize
    #a.qual= r.qual
    #a.tags = r.tags

    #if len(items) == 6:
    #    a.rname = items[1]
    #    a.pos += int( items[3] )
    #    a.mpos += int( items[3] )
    #out.write( a )

out.close()
f.close()
