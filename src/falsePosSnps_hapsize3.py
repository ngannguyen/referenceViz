#!/usr/bin/env python

"""
nknguyen soe ucsc edu
May 22 2012
This script examines potential false positive SNPs in the n->1 problem, where duplications exist in the population but not the reference, and reads of the absent copy falsely mapped to the copy the reference has.

Particularly, we will use trio data and linkage to detect if a SNP is more likely to be false or not
Input: a/ file that contains regions of interest
       b/ phased trio SNPs calls (from 1k GP), which should include SNPs calls for the father, mother and child
Output: 
"""

import os, re, sys


class TrioSnp():
    def __init__(self, line, ref, alts):
        alleleStr = line.split(':')[0]
        
        self.phased = False
        items = alleleStr.split('/')
        if re.search("\|", alleleStr):
            self.phased = True
            items = alleleStr.split('|')
        if len(items) != 2:
            raise ValueError("SNP does not have two alleles: %s" % line)
        
        allAlleles = [ref]
        allAlleles.extend(alts)
        self.alleles = [ allAlleles[int(i)] for i in items ]

class Snp():
    def __init__(self, line):
        line = line.strip()
        items = line.split('\t')
        self.chr = items[0].lstrip('chr')
        self.pos = int(items[1])
        self.id = items[2]
        self.ref = items[3]
        self.alts = items[4].split(',')
        self.qual = items[5]
        self.filter = items[6]
        self.father = TrioSnp(items[9], self.ref, self.alts)
        self.mother = TrioSnp(items[10], self.ref, self.alts)
        self.child = TrioSnp(items[11], self.ref, self.alts)

def getPhasedSnps(pos2snp):
    pos2phasedSnp = {}
    for pos, snp in pos2snp.iteritems():
        if snp.father.phased and snp.mother.phased and snp.child.phased:
            pos2phasedSnp[pos] = snp
    return pos2phasedSnp

def checkTrioSnps(chr2pos, chr2pos2snp):
    #chr2pos is a list of suspicious positions
    #chr2pos2snp is the list of 1k genome snps of the trios (father, mother and child)
    #for each suspicious position, check to see if there is a phased snp called by 1kGP
    #If yes, look for the closest left phased snp, and closest right phased snp
    #and check to see if the child haplotype requires any recombination, if yes and the 
    #distance between left, current, and right SNPs are very close, it is probably a spurious SNP
    
    numRecoms = 0 #number of cases where recombination is needed to explain 
    numExams = 0 #number of snp tested
    range = 1000
    for chr in chr2pos:
        if chr not in chr2pos2snp:
            continue
        pos2snp = chr2pos2snp[chr]
        pos2phasedSnps = getPhasedSnps(pos2snp) #phased SNPs = SNPs where both parents & child are phased
        sortedPos = sorted(pos2phasedSnps.keys())
        pos2index = {}
        for i, p in enumerate(sortedPos):
            pos2index[p] = i

        for pos in chr2pos[chr].keys():
            if pos not in pos2phasedSnps:
                continue
            snp = pos2snp[pos]
            if snp.child.alleles[0] == snp.ref and snp.child.alleles[1] == snp.ref: #no SNP in the child
                continue

            i = pos2index[pos]
            if i == 0 or i == len(pos2index):
                continue

            lefti = i - 1
            leftpos = sortedPos[lefti]
            #while leftpos in chr2pos[chr].keys() and pos - leftpos <= 1000 and lefti > 0:
            #    lefti -= 1
            #    leftpos = sortedPos[lefti]
            #if leftpos in chr2pos[chr].keys() or pos - leftpos > 1000:
            #    continue

            righti = i + 1
            rightpos = sortedPos[righti]
            #while rightpos in chr2pos[chr].keys() and rightpos - pos <= 1000 and righti < len(sortedPos) - 1:
            #    righti += 1
            #    rightpos = sortedPos[righti]
            #if rightpos in chr2pos[chr].keys() or rightpos - pos > 1000:
            #    continue
            
            left = pos2phasedSnps[ leftpos ]
            right = pos2phasedSnps[ rightpos ]
            numExams += 1
            
            childhaps = [] 
            for i in [0, 1]:
                hapi = [left.child.alleles[i], snp.child.alleles[i], right.child.alleles[i]]
                childhaps.append(hapi)
            
            motherhaps = [] 
            for i in [0, 1]:
                hapi = [left.mother.alleles[i], snp.mother.alleles[i], right.mother.alleles[i]]
                motherhaps.append(hapi)
            
            fatherhaps = [] 
            for i in [0, 1]:
                hapi = [left.father.alleles[i], snp.father.alleles[i], right.father.alleles[i]]
                fatherhaps.append(hapi)
            
            print left.pos, snp.pos, right.pos
            print motherhaps
            print fatherhaps
            print childhaps

            ok = False
            if (childhaps[0] in motherhaps and childhaps[1] in fatherhaps) or (childhaps[0] in fatherhaps and childhaps[1] in motherhaps):
                ok = True
            if not ok:
                numRecoms += 1
    sys.stdout.write("Number of SNPs examed: %d\n" %numExams)
    sys.stdout.write("Potential spurious SNPs: %d\t%f%%\n" %(numRecoms, numRecoms*100.0/numExams))

def readPosFile(file):
    f = open(file, 'r')
    chr2pos = {}
    for line in f:
        items = line.split('\t')
        chritems = items[0].split('.')
        #hg18.chr6.170899992.28585733.4970599.1
        chr = chritems[1].lstrip('chr')
        offset = int(chritems[3])
        pos = int(items[1]) + offset
        if chr not in chr2pos:
            chr2pos[chr] = {pos: 1}
        else:
            chr2pos[chr][pos] = 1
    f.close()
    return chr2pos

def readVcfFile(file):
    chr2pos2snp = {}
    f = open(file, 'r')
    for line in f:
        if line[0] == '#':
            continue
        snp = Snp(line)
        if snp.chr not in chr2pos2snp:
            chr2pos2snp[snp.chr] = { snp.pos: snp }
        else:
            chr2pos2snp[snp.chr][snp.pos] = snp
    f.close()
    return chr2pos2snp 

def main():
    vcffile = sys.argv[1]
    posfile = sys.argv[2]
    chr2pos = readPosFile(posfile)
    chr2pos2snp = readVcfFile(vcffile)
    checkTrioSnps(chr2pos, chr2pos2snp)


if __name__ == '__main__':
    main()


