#!/usr/bin/env python

"""
Compare Snps between each sample and the reference with the set of Snps in dbSnps
"""

import os, sys, re
#import sqlite3
from optparse import OptionParser
import xml.etree.ElementTree as ET
from sonLib.bioio import system

class SnpSite():
    def __init__(self, line, sampleName, referenceName):
        items = line.strip().split()
        self.name = items[0]
        self.start = int( items[1] )
        self.allele = items[2]
        self.ref = items[3]
        self.refstart = int( items[4] ) 
        refitems = self.ref.split('.')
        if len(refitems) <=2:
            sys.stderr.write('Reference sequence must provide chromosome info: i.e: species.chr\n')
            sys.exit(1)
        self.refchrom = refitems[1]

        #Convert coordinate if necessary:
        if len(refitems) == 6:
            refstrand = refitems[5]
            reflen = int(refitems[2])
            offset = int(refitems[3])
            fraglen = int(refitems[4])
            if refstrand == '1':#Positive strand:
                self.refstart = self.refstart + offset
            else:
                self.refstart = self.refstart + reflen - (offset + fraglen)
        self.refallele = items[5]
        self.sampleName = sampleName
        self.referenceName = referenceName
    
    def __cmp__(self, other):
        chr = int( self.refchrom.lstrip('chr') )
        otherchr = int( other.refchrom.lstrip('chr') )

        if chr < otherchr:
            return -1
        elif chr > otherchr:
            return 1
        else:
            if self.refstart < other.refstart:
                return -1
            elif self.refstart == other.refstart:
                return 0
            else:
                return 1
    #def getDbFormat(self):
    #    return (self.name, self.start, self.allele, self.ref, self.refstart, self.refallele, self.sampleName, self.referenceName)

class Snp():
    def __init__(self, line):
        items = line.split('\t')
        if len(items) < 25:
            sys.stderr.write("Each snp record must have at least 25 fields\n")
            sys.stderr.write("%s\n" %line)
        self.chrom = items[0]
        self.chromStart = int(items[1])
        self.chromEnd = int(items[2])
        self.name = items[3]
        self.score = int(items[4])
        self.strand = items[5]
        self.refNCBI = items[6]
        self.refUCSC = items[7]
        self.observed = items[8]
        self.molType = items[9]
        self.type = items[10]
        #self.valid = items[11]
        self.func = items[14]
        self.locType = items[15]
        self.exceptions = items[17]
        self.alleleFreqCount = int(items[20])
        self.alleles = items[21].rstrip(',')
        self.alleleFreqs = items[23].rstrip(',')
        #self.alleles = items[21].rstrip(',').split(',')
        #self.alleleFreqs = [ float(f) for f in items[23].rstrip(',').split(',') ]

    def __cmp__(self, other):
        chr = int( self.chrom.lstrip('chr') )
        otherchr = int( other.chrom.lstrip('chr') )

        if chr < otherchr:
            return -1
        elif chr > otherchr:
            return 1
        else:
            if self.chromStart < other.chromStart:
                return -1
            elif self.chromStart == other.chromStart:
                return 0
            else:
                return 1
    #def reverse(self):#convert to the opposite strand:
    #    chrSizes = {'chr6':171115067}
    #    chrsize = chrSizes[self.chrom]
    #    if self.strand == '-':
    #        self.strand = '+'
    #    else:
    #        self.strand = '-'
    #    self.chromStart = 
    
    #def getDbFormat(self):
    #    return (self.chrom, self.chromStart, self.chromEnd, self.name, self.score,\
    #            self.strand, self.refNCBI, self.refUCSC, self.observed, self.molType,\
    #            self.type, self.func, self.locType, self.exceptions, self.alleleFreqCount, self.alleles, self.alleleFreqs)
    
def readDbSnps(file):
    f = open(file, 'r')
    snps = []
    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        if snp.type == 'single':
            snps.append( snp )

    return snps

def readRefSnps(file, filteredSamples):
    xmltree = ET.parse( file )
    root = xmltree.getroot()
    snps = {}
    samples = []
    for sample in root.findall( 'statsForSample' ):
        name = sample.attrib['sampleName']
        if name == 'ROOT' or name == '' or name in filteredSamples:
            continue
        samples.append(name)
        snps[name] = []
        ref = sample.attrib['referenceName']
        #totalSnps = int( sample.attrib['totalErrors'] )
        sites = sample.text.strip().split('\n')
        for site in sites:#each snp detected by cactus ref, check to see if there is one in dbsnp
            if site != '':
               snp = SnpSite(site, name, ref)
               snps[name].append( snp )
    return snps, samples

def getStats(dbsnps, refsnps, samples):
    dbsnps.sort()
    #print [ s.chromStart for s in dbsnps]

    totalSnps = len(dbsnps)
    sys.stdout.write("Sample\tTotalCalledSnps\tTP\tPercentageTP\tFP\tPercentageFP\n")
    for sample in samples:
        #sys.stderr.write("Sample %s, sorting...\n" %(sample))
        snps = sorted( refsnps[sample] )
        #sys.stderr.write("Done sorting, len %d\n" %(len(snps)))
        #print [s.refstart for s in snps]

        refTotal = len(snps)
        tp = 0
        currindex = 0

        for s in snps:
            #flag = False
            #print "CURRENT INDEX: %d" %currindex
            for i in xrange(currindex,totalSnps):
                dbs = dbsnps[i].chromStart
                currindex = i
                if s.refstart > dbs:
                    #print "\tnot there yet: %d < %d" %(s.refstart, dbs)
                    continue
                elif s.refstart == dbs:
                    #print "\tfound it! %d, %d" %(s.refstart, dbs)
                    tp +=1
                    #flag = True
                    break
                else:
                    #print "\tGone too far, stop: %d, %d" %(s.refstart, dbs)
                    break
            #if flag == False:
            #    sys.stderr.write("%s\t%d\n" %(sample, s))

        fp = refTotal - tp
        fn = totalSnps - tp

        sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\n" %(sample, refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))

def initOptions( parser ):
    #parser.add_option('-o', '--outfile', dest='outfile', default='', help='Output file. Default is stdout' )
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def checkOptions( args, options, parser ):
    if len(args) < 2:
        parser.error('Please provide two input files: cactusRefSnpFile and dbSnpFile\n')
    #options.snpFile = args[0]
    if not os.path.exists(args[0]):
        parser.error('File %s does not exist\n' %args[0])
    #options.dbSnpFile = args[1]
    if not os.path.exists(args[1]):
        parser.error('File %s does not exist\n' %args[1])
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split('-')
    else:
        options.filteredSamples = []

def main():
    usage = ('Usage: %prog [options] snpStats_*.xml snp134dump.txt')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    
    dbsnps = readDbSnps( args[1] )
    refsnps,samples = readRefSnps( args[0], options.filteredSamples )
    getStats( dbsnps, refsnps, samples )

    #Delete dbfile, refdbfile ...

if __name__ == '__main__':
    main()

