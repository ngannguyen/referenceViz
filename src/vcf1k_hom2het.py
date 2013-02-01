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

import os, re, sys, random
import numpy as np


class Sample():
    def __init__(self, line, ref, alts):
        alleleStr = line.split(':')[0]
        self.ref = ref
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
        self.het = True
        self.homref = False
        self.homalt = False
        if self.alleles[0] == self.alleles[1]:
            self.het = False
            if self.alleles[0] == ref:
                self.homref = True
            else:
                self.homalt = True

class Snp():
    def __init__(self, line):
        line = line.strip()
        items = line.split('\t')
        self.chr = items[0].lstrip('chr')
        self.pos = int(items[1])
        self.id = items[2]
        self.ref = items[3]
        self.alts = items[4].split(',')
        
        isindel = False
        if len(self.ref) > 1:
            isindel = True
        else:
            for a in self.alts:
                if len(a) > 1:
                    isindel = True
                    break
        self.isindel = isindel

        self.qual = items[5]
        self.filter = items[6]
        samples = []
        for item in items[9:]:
            if item == '' or len( item.split(':') ) < 2:
                continue
            samples.append( Sample(item, self.ref, self.alts) )
        self.hetfreq = getHetFreq(samples)
        self.homfreq = 1 - self.hetfreq
        self.homreffreq = getHomRef(samples)
        self.homaltfreq = 1 - self.hetfreq - self.homreffreq

def getHomRef(samples):
    total = len(samples)
    numhomref = 0
    for s in samples:
        if s.homref:
            numhomref += 1
    return numhomref*1.0/total

def getHetFreq(samples):
    total = len(samples)
    numhets = 0
    for s in samples:
        if s.het:
            numhets += 1
    return numhets*1.0/total

def homhetstatsSampling(pos2snp, size):
    snps = random.sample( pos2snp.values(), size )
    total = len(snps)
    het = 0
    homref = 0
    homalt = 0
    
    for snp in snps:
        het += snp.hetfreq
        homref += snp.homreffreq
        homalt += snp.homaltfreq
    return homref/total, het/total, homalt/total

def sampling(pos2snp, numsamplings, size):
    hets = []
    homrefs = []
    homalts = []
    for i in xrange(numsamplings):
        homref, het, homalt = homhetstatsSampling( pos2snp,  size )
        hets.append(het)
        homrefs.append(homref)
        homalts.append(homalt)
    return np.mean(homrefs), np.std(homrefs), np.mean(hets), np.std(hets), np.mean(homalts), np.std(homalts)

def homhetstats(chr2pos, chr2pos2snp):
    hets = []
    homrefs = []
    homalts = []
    selectedhets = []
    selectedhomrefs = []
    selectedhomalts = []

    for chr in chr2pos2snp:
        for pos, snp in chr2pos2snp[chr].iteritems():
            hets.append( snp.hetfreq )
            homrefs.append( snp.homreffreq )
            homalts.append( snp.homaltfreq )
            if chr in chr2pos and pos in chr2pos[chr]:
                selectedhets.append( snp.hetfreq )
                selectedhomrefs.append( snp.homreffreq )
                selectedhomalts.append( snp.homaltfreq )
    hetavr = np.mean( hets )
    hetstd = np.std( hets )
    homrefavr = np.mean( homrefs )
    homrefstd = np.std( homrefs )
    homaltavr = np.mean( homalts )
    homaltstd = np.std( homalts )
    
    selectedhetavr = np.mean( selectedhets )
    selectedhetstd = np.std( selectedhets )
    selectedhomrefavr = np.mean( selectedhomrefs )
    selectedhomrefstd = np.std( selectedhomrefs )
    selectedhomaltavr = np.mean( selectedhomalts )
    selectedhomaltstd = np.std( selectedhomalts )

    print "Category\tAvr AA\tAvr AB\tAvr BB\tStd AA\tStd AB\tStd BB"
    print "Overall\t%f\t%f\t%f\t\t%f\t%f\t%f" %(homrefavr, hetavr, homaltavr, homrefstd, hetstd, homaltstd)
    print "Selected\t%f\t%f\t%f\t\t%f\t%f\t%f" %(selectedhomrefavr, selectedhetavr, selectedhomaltavr, selectedhomrefstd, selectedhetstd, selectedhomaltstd)
    print "Selected/Overall\t%f\t%f\t%f" %(selectedhomrefavr/homrefavr, selectedhetavr/hetavr, selectedhomaltavr/homaltavr)

    numsamplings = 100
    #size = 4
    size = 1000
    print "\nSampling, size %d, %d times" %(size, numsamplings)
    samplinghomrefavr, samplinghomrefstd, samplinghetavr, samplinghetstd, samplinghomaltavr, samplinghomaltstd = sampling(chr2pos2snp['6'], numsamplings, size)
    print "Sampling\t%f\t%f\t%f\t\t%f\t%f\t%f" %(samplinghomrefavr, samplinghetavr, samplinghomaltavr, samplinghomrefstd, samplinghetstd, samplinghomaltstd)

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
        sys.stderr.write('Done %s\n' %snp.pos)
        if snp.filter != 'PASS' or snp.isindel:
            continue
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
    #checkTrioSnps(chr2pos, chr2pos2snp)
    #checkUnphasedSnps(chr2pos, chr2pos2snp)
    #countHaps(chr2pos, chr2pos2snp)
    homhetstats(chr2pos, chr2pos2snp)


if __name__ == '__main__':
    main()


