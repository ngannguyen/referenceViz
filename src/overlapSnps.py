#!/usr/bin/env python

'''
Oct 20 2011: nknguyen soe ucsc edu
Take an intersect of two input snps file
Each input snp file has the format:
Sample\tchrom\tstart\tend\taltAllele\t...
'''

import re, os, sys

class Snp():
    def __init__(self, line):
        self.desc = line
        items = line.strip().split('\t')
        if len(items) < 5:
            sys.stderr.write("Wrong input format, required at least 5 fields, only got %d. Line: %s\n" %(len(items), line))
            sys.exit(1)
        self.sample = items[0]
        self.chrom = items[1].lstrip('chrom').lstrip('chr')
        self.start = int(items[2])
        self.end = int(items[3])
        self.allele = items[4].upper()

def readfile( file ):
    snps = {}
    f = open(file, 'r')
    for line in f:
        snp = Snp(line)
        if snp.sample not in snps:
            snps[snp.sample] = [snp]
        else:
            snps[snp.sample].append(snp)
    return snps
    f.close()

def overlap( snps1, snps2, file):
    f = open(file, 'w')
    sys.stdout.write("Sample\tcactusSnps\tmpileupSnps\toverlap\tPercentage\n")
    for sample in snps1:
        if sample not in snps2:
            continue
        list1 = snps1[sample]
        list2 = snps2[sample]
        overlapCount = 0

        i = 0
        j = 0

        for s1 in list1:
            i = j
            flag = False
            if i < len(list2):
                s2 = list2[i]
                while s2.start <= s1.start:
                    if s2.start == s1.start and s2.end == s1.end and s1.allele == s2.allele:
                        flag = True
                        break
                    if s2.start < s1.start:
                        j += 1
                    i +=1
                    if i == len(list2):
                        break
                    s2 = list2[i]
            if flag:
                overlapCount += 1
            else:
                f.write('%s' % s1.desc)
        sys.stdout.write("%s\t%d\t%d\t%d\t%.2f%%\n" %(sample, len(list1), len(list2), overlapCount, overlapCount*100.0/len(list1)))
    f.close()

snps1 = readfile( sys.argv[1] )
snps2 = readfile( sys.argv[2] )

overlap(snps1, snps2, "%s-not-%s" %(os.path.basename(sys.argv[1]).split('.')[0], os.path.basename(sys.argv[2]).split('.')[0]))

