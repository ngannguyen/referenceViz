#!/usr/bin/env python

"""
Compare Snps between each sample and the reference with the set of Snps in dbSnps
"""

import os, sys, re
#import sqlite3
from optparse import OptionParser
import xml.etree.ElementTree as ET
from sonLib.bioio import system

class PgSnp():
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 8:
            sys.stderr.write("pgSnp table, line '%s' does not have enough fields. 8 is expected, only got %d\n" %(line, len(items)))
            sys.exit(1)
        self.chrom = items[0]
        self.chromStart = int(items[1])
        self.chromEnd = int(items[2])
        self.name = items[3]
        self.alleles = [ a.lower() for a in self.name.split('/') ]
        self.sample = items[7]
        
        #Check if it's snp or not:
        self.isSnp = False
        for o in self.alleles:
            if len(o) == 1 and o != '-':
                self.isSnp = True
                break
        
        #self.alleleCount = int(items[4])
        #self.alleleFreq = [ int(f) for f in items[5].split(',') ]
        #if len(self.alleleFreq) != self.alleleCount:
        #    sys.stderr.write("")

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

class SnpSite():
    def __init__(self, line, sampleName, referenceName):
        items = line.strip().split()
        self.name = items[0]
        self.start = int( items[1] )
        self.allele = items[2].lower()
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
        self.refallele = items[5].lower()
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
        self.refNCBI = items[6].lower()
        self.refUCSC = items[7].lower()
        self.observed = items[8].lower().split('/')
        if self.strand == '-':
        #    self.refNCBI = reverse(self.refNCBI)
        #    self.refUCSC = reverse(self.refUCSC)
            self.observed = [ reverse(a) for a in self.observed] 
        #Check if indel is snp:
        self.isSnp = False
        for o in self.observed:
            if len(o) == self.chromEnd - self.chromStart and o!= '-' and o != self.refUCSC:
                self.isSnp = True
                break
        
        self.molType = items[9]
        self.type = items[10]
        #self.valid = items[11]
        self.func = items[14]
        self.locType = items[15]
        self.exceptions = items[17]
        self.alleleFreqCount = int(items[20])
        #self.alleles = items[21].rstrip(',')
        #self.alleleFreqs = items[23].rstrip(',')
        self.alleles = items[21].lower().rstrip(',').split(',')
        #for allele in self.alleles:
        #    if len(allele) < self.chromEnd - self.chromStart:
        #        sys.stderr.write("dbsnp line '%s', allele length is different from chromEnd - chromStart\n" %(line))
        #        sys.exit(1)
        if self.strand == '-':
            d = {'a':'t', 't':'a', 'c':'g', 'g':'c'}
            newalleles = []
            for a in self.alleles:
                if a in d:
                    newalleles.append( d[a] )
                else:
                    newalleles.append( a )
            self.alleles = newalleles

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

def reverse(s):
    d = {'a':'t', 't':'a', 'c':'g', 'g':'c'}
    if s in d:
        return d[s]
    else:
        return s

def readDbSnps(file):
    f = open(file, 'r')
    snps = []

    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        #if snp.type == 'single' or snp.type == 'mnp':
        if snp.isSnp:
            snps.append( snp )
    sys.stdout.write("TotalSnps\t%d\n" % (len(snps)))

    return snps

def readPgSnp(file):
    f = open(file, 'r')
    sample2snps = {}
    for line in f:
        if re.search(line, 'chromStart'):
            continue
        snp = PgSnp(line)
        if not snp.isSnp:
            continue
        if snp.sample not in sample2snps:
            sample2snps[snp.sample] = [snp]
        else:
            sample2snps[snp.sample].append(snp)
    f.close()
    return sample2snps


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

def checkAlleles(allele, refallele, snp):
    if refallele != snp.refNCBI and refallele != snp.refUCSC:
        return False
    elif allele not in snp.observed:
        return False
    return True

def checkPgAlleles(allele, pgsnp):
    if allele not in pgsnp.alleles:
        return False
    return True
   
def calcOverlapSnps(snps, reportedSnps, isPgSnp):#reportedSnps can be dbSnps or pgSnps
    refTotal = len(snps)
    totalSnps = len(reportedSnps)
    tpPos = 0 #pass position check (coordinate check) 
    tp = 0 #pass both coordinate and allele check
    currindex = 0

    for s in snps:
        flag = False
        for i in xrange(currindex,totalSnps):
            dbs = reportedSnps[i].chromStart
            dbe = reportedSnps[i].chromEnd

            if flag == False:
                currindex = i
            if dbe <= s.refstart :
                #print "\tnot there yet: %d < %d" %(s.refstart, dbs)
                continue
            elif dbs <= s.refstart : # < dbe
                flag = True
                #print "\tfound it! %d, %d" %(s.refstart, dbs)
                #tpPos += 1
                if (isPgSnp and checkPgAlleles(s.allele, reportedSnps[i]) ) or (not isPgSnp and checkAlleles(s.allele, s.refallele, reportedSnps[i])):
                #if checkAlleles(s.allele, s.refallele, reportedSnps[i].alleles):
                    tp += 1
                    break
                #else:
                #    sys.stderr.write("%s\t%d\t%d\tourref: %s\tourAllele: %s\tncbi: %s\tucsc: %s\tobserved: %s\n" %(s.sampleName, s.refstart, dbs, s.refallele, s.allele, reportedSnps[i].refNCBI, reportedSnps[i].refUCSC, ','.join(reportedSnps[i].observed) ))
            else: #s.refstart < dbs
                #print "\tGone too far, stop: %d, %d" %(s.refstart, dbs)
                break
        if flag == True: #pass position check
            tpPos += 1

    return tp, tpPos

def getStats(dbsnps, refsnps, samples, sample2snps):
    dbsnps.sort()
    #print [ s.chromStart for s in dbsnps]

    totalSnps = len(dbsnps)
    sys.stdout.write("Sample\tTotalCalledSnps\ttpPos\tPercentageTpPos\tTP\tPercentageTP\tFP\tPercentageFP\tsampleSnps\tsampleTpPos\tPercentageSampleTpPos\tsampleTP\tPercentageSampleTP\tsampleFN\tPercentageSampleFN\n")
    for sample in samples:
        #sys.stderr.write("Sample %s, sorting...\n" %(sample))
        snps = sorted( refsnps[sample] )
        #sys.stderr.write("Done sorting, len %d\n" %(len(snps)))
        
        #Check against dbSnps
        tp, tpPos = calcOverlapSnps(snps, dbsnps, False)
        refTotal = len(snps)
        fp = refTotal - tp 
        if refTotal > 0:
            #sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
            sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tpPos, 100.0*tpPos/refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
        else:
            sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tpPos, 0.00, tp, 0.00, fp, 0.00))

        #check against snps that were reported specifically for the sample
        if sample in sample2snps:
            pgsnps = sample2snps[sample]
            pgsnps.sort()
            pgTp, pgTpPos = calcOverlapSnps(snps, pgsnps, True)
            pgTotal = len(pgsnps)
            fn = pgTotal - pgTp
            if refTotal > 0 and pgTotal > 0:
                #sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
                sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 100.0*pgTpPos/refTotal, pgTp, 100.0*pgTp/refTotal, fn, 100.0*fn/pgTotal))
            else:
                sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 0.00, tp, 0.00, fn, 0.00))

        #check against snps that were reported specifically for the sample

        sys.stdout.write("\n")


def initOptions( parser ):
    #parser.add_option('-o', '--outfile', dest='outfile', default='', help='Output file. Default is stdout' )
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    parser.add_option('--pgSnp', dest='pgSnp', help="pgSnp file (snps found in each sample)")

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
    if options.pgSnp and not os.path.exists(options.pgSnp):
        parser.error("pgSnp file %s does not exists." %(options.pgSnp))

def main():
    usage = ('Usage: %prog [options] snpStats_*.xml snp134dump.txt')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    
    dbsnps = readDbSnps( args[1] )
    sample2snps = {}
    if options.pgSnp:
        sample2snps = readPgSnp(options.pgSnp)
    refsnps,samples = readRefSnps( args[0], options.filteredSamples )
    getStats( dbsnps, refsnps, samples, sample2snps )

    #Delete dbfile, refdbfile ...

if __name__ == '__main__':
    main()

