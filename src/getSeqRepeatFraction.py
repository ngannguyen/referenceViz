#!/usr/bin/env python

"""
Get fraction of input sequences that are masked by repeatMasker using faSize
"""
import os, sys, re, copy
from sonLib.bioio import system

indir = sys.argv[1]
seqStr = "NA12878  NA12892  NA19238  NA19239  NA19240  apd  cox  dbb  hg19  mann  mcf  nigerian  panTro3  qbl  ssto  venter  yanhuang"
seqs = seqStr.split()

#Parse faSizeOut to get repeat fraction:
def readFaSizeOut(file):
    f = open(file, 'r')
    line = f.readline()
    #4809387 bases (0 N's 4809387 real 2289302 upper 2520085 lower) in 128 sequences in 1 files
    items = line.split()
    if len(items) < 16:
        sys.stderr.write("faSizeOut %s has wrong format\n" %file)
        sys.exit(1)
    total = int(items[0])
    repeat = int(items[8])
    repeatPercentage = 0.0
    if total > 0:
        repeatPercentage = repeat*100.0/total
    return total, repeat, repeatPercentage

#Call faSize:
sys.stdout.write("Sample\tTotal\tRepeatBases\tRepeatPercentage\n")
for seq in seqs:
    infile = os.path.join(indir, seq)
    if not os.path.exists(infile):
        continue

    faSizeOut = "%s-faSizeOut.txt" %seq
    system("faSize %s > %s" %(infile, faSizeOut))
    total, repeat, repeatPc = readFaSizeOut( faSizeOut )
    sys.stdout.write("%s\t%d\t%d\t%.2f\n" %(seq, total, repeat, repeatPc))
    system("rm %s" %faSizeOut)

