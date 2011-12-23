#!/usr/bin/env python

import os, sys, re
import pysam

file1 = sys.argv[1]
file2 = sys.argv[2]
outdir = sys.argv[3]
f1 = pysam.Samfile(file1, "rb") #Input bam from stdin
f2 = pysam.Samfile(file2, "rb") #Input bam from stdin

def decodeFlag(flag):
    flagDict = {"isPaired": 0, "properlyPaired":1, "unmapped":2, \
                "unmappedMate":3, "strand":4, "mateStrand":5, "firstRead":6,\
                "secondRead":7, "notPrimaryAlign":8, "failedQual":9,\
                "duplicate":10}

    for key in flagDict:
        if 2**flagDict[key] & flag == 0:
            flagDict[key] = False
        else:
            flagDict[key] = True
    return flagDict

#===========
def readBam( f ):
    uniqs = []
    multi = []
    for r in f:
        flagDict = decodeFlag( r.flag )
        if flagDict["unmapped"]:
            continue
        name = r.qname
        if flagDict["isPaired"]:
            if flagDict["firstRead"]:
                name = name + "_1"
            else:
                name = name + "_2"
        r.qname = name
        if r.tags != None and 'XA' in [ tag for (tag, vals) in r.tags ]:#not uniq
            #multi.append(name)
            multi.append(r)
        else:
            #uniqs.append(name)
            uniqs.append(r)
    return uniqs, multi

uniqs1 , multi1 = readBam( f1 )
sys.stderr.write("File %s has %d uniquelyMapped reads, and %d multi-mapped reads\n" %(file1, len(uniqs1), len(multi1)))
uniqs2 , multi2 = readBam( f2 )
sys.stderr.write("File %s has %d uniquelyMapped reads, and %d multi-mapped reads\n" %(file2, len(uniqs2), len(multi2)))

uniqs1 = sorted(uniqs1, key=lambda read:read.qname )
sys.stderr.write("Sorted uniqs1\n")
multi1 = sorted(multi1, key=lambda read:read.qname)
sys.stderr.write("Sorted multi1\n")
uniqs2 = sorted(uniqs2, key=lambda read:read.qname)
sys.stderr.write("Sorted uniqs2\n")
multi2 = sorted(multi2, key=lambda read:read.qname)
sys.stderr.write("Sorted multi2\n")

def seperateOverlap( list1, list2 ):
    overlaps = []
    u1 = []
    u2 = []
    i = 0
    j = 0
    while i < len( list1 ):
        r1 = list1[i]
        if j >= len(list2):
            u1.append(r1)
            i += 1
        else:
            r2 = list2[j]
            if r2.qname > r1.qname:
                u1.append(r1)
                i += 1
            else:
                while r2.qname <= r1.qname:
                    if r2.qname == r1.qname:
                        overlaps.append(r2)
                        j += 1
                        i += 1
                        break
                    else:
                        u2.append(r2)
                        j+= 1
                    if j >= len(list2):
                        break
                    r2 = list2[j]
    return overlaps, u1, u2

uniqs_overlaps, u1, u2 = seperateOverlap( uniqs1, uniqs2 )
sys.stderr.write("%d uniqs_overlap, %d u1, %d u2\n" %(len(uniqs_overlaps), len(u1), len(u2)))
multi_overlaps, m1, m2 = seperateOverlap( multi1, multi2 )
sys.stderr.write("%d multi_overlap, %d m1, %d m2\n" %(len(multi_overlaps), len(m1), len(m2)))

#sys.stdout.write("Uniquely mapped in %s, Multimapped in %s\n" %(file1, file2))

file1name = os.path.basename(file1).rstrip('.bam')
file2name = os.path.basename(file2).rstrip('.bam')
out1m = pysam.Samfile(os.path.join(outdir, "%sM-%sU_M.bam" %(file2name, file1name)), "wb", template=f2)
out1u = pysam.Samfile(os.path.join(outdir, "%sM-%sU_U.bam" %(file2name, file1name)), "wb", template=f1)
for r2 in m2:
    for r1 in u1: 
    #if r.qname in [r1.qname for r1 in u1]:
        if r2.qname == r1.qname:
            out1m.write(r2)
            out1u.write(r1)
            break

out2m = pysam.Samfile(os.path.join(outdir, "%sM-%sU_M.bam" %(file1name, file2name)), "wb", template=f1)
out2u = pysam.Samfile(os.path.join(outdir, "%sM-%sU_U.bam" %(file1name, file2name)), "wb", template=f2)
for r1 in m1:
    for r2 in u2:
        if r1.qname == r2.qname:
            out2m.write(r1)
            out2u.write(r2)
            break
#sys.stdout.write("Multimapped in %s, uniquely mapped in %s\n" %(file1, file2))
#for r in m1:
#    if r in u2:
#        sys.stdout.write("%s\n" %r)

f1.close()
f2.close()
out1m.close()
out1u.close()
out2m.close()
out2u.close()

