#!/usr/bin/env python

"""
Compare Snps between each sample and the reference with the set of Snps in dbSnps
"""

import os, sys, re, copy
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
    def __init__(self, line, sampleName, referenceName, longnames):
        items = line.strip().split()
        self.name = items[0]
        self.start = int( items[1] )
        self.allele = items[2].lower()
        self.ref = items[3]
        if longnames and self.ref in longnames:
            self.ref = longnames[self.ref]

        self.refstart = int( items[4] ) 
        refitems = self.ref.split('.')
        if len(refitems) <=2:
            self.refchrom = self.ref
            #sys.stderr.write('Reference sequence must provide chromosome info: i.e: species.chr\n')
            #sys.exit(1)
        else:
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
        #chr = int( self.refchrom.lstrip('chr') )
        #otherchr = int( other.refchrom.lstrip('chr') )
        chr = self.refchrom.lstrip('chr')
        otherchr = other.refchrom.lstrip('chr')

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
        self.observed.sort()
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
        self.alleles = items[21].lower().rstrip(',').split(',')
        self.alleleFreqs = {}
        if items[23] != "":
            alleleFreqs = [ float(freq) for freq in items[23].rstrip(',').split(',') ]
            for i, allele in enumerate(self.alleles):
                self.alleleFreqs[ allele ] = alleleFreqs[i]
        #self.alleleFreqCount = int(items[20])
        #self.alleles = items[21].rstrip(',')
        #self.alleleFreqs = items[23].rstrip(',')
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
        #chr = int( self.chrom.lstrip('chr') )
        #otherchr = int( other.chrom.lstrip('chr') )
        chr = self.chrom.lstrip('chr')
        otherchr = other.chrom.lstrip('chr')

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

class Filter():
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 3:
            sys.stderr.write("Filter region has wrong format: %s\n" %line)
        self.chrom = items[0]
        self.chromStart = int(items[1])
        self.chromEnd = int(items[2])

def splitDbSnp(snp):
    snps = []
    #if len(snp.observed) == 0 or len(snp.refNCBI) != len(snp.observed[0]) or len(snp.refUCSC) != len(snp.observed[0]):
    if len(snp.observed) == 0:
        return snps
    snpLen = len(snp.observed[0])
    for allele in snp.observed:
        if len(allele) != snpLen:
            return snps
    #if len(snp.observed) == 0 or len(snp.refNCBI) != len(snp.observed[0]) or len(snp.refUCSC) != len(snp.observed[0]):
    #    #sys.stderr.write("%d, Len refNCBI %s %d, len refUCSC %s %d, len observed[0]: %s %d\n" %( snp.chromStart, snp.refNCBI, len(snp.refNCBI), snp.refUCSC, len(snp.refUCSC), snp.observed[0], len(snp.observed[0])))
    #    return snps
    #sys.stderr.write("Splitting snps...\n")
    for i in xrange(snpLen):
        s = copy.copy( snp )
        if len(snp.refNCBI) > i:
            s.refNCBI = snp.refNCBI[i]
        else:
            s.refNCBI = ''
        if len(snp.refUCSC) > i:
            s.refUCSC = snp.refUCSC[i]
        else:
            s.refUCSC = ''
        s.chromStart = snp.chromStart + i
        s.observed = [ a[i] for a in snp.observed ]
        s.isSnp = True
        s.type = 'single'
        snps.append(s)
    return snps

def reverse(s):
    d = {'a':'t', 't':'a', 'c':'g', 'g':'c'}
    if s in d:
        return d[s]
    else:
        return s

def isInRange(coord, start, end):
    if start != None and coord < start:
        return False
    if end != None and coord > end:
        return False
    return True

def isInFilter(coord, filter):
    for f in filter:
        if f.chromStart <= coord and coord < f.chromEnd:
            return True
    return False

def readFilter(file):
    f = open(file, 'r')
    regions = []
    for line in f:
        if re.search('tart', line):
            continue
        regions.append( Filter(line) )
    f.close()
    return regions

def readSeqNameFile(file):
    names = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split("\t")
        if len(items) < 2:
            continue
        assert items[1] not in names
        names[items[1]] = items[0]
    f.close()
    return names

def readDbSnps(file, start, end, filter):
    f = open(file, 'r')
    snps = []
    nummnp = 0
    i = 0
    j = 0

    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        #if snp.type == 'single' or snp.type == 'mnp':
        inrange = isInRange(snp.chromStart, start, end)
        #infilter = isInFilter(snp.chromStart, filter)
        
        infilter = False
        if len(filter) > i:
            i = j
            filterReg = filter[i]
            while filterReg.chromStart <= snp.chromStart:
                if filterReg.chromStart <= snp.chromStart and snp.chromStart < filterReg.chromEnd:
                    infilter = True
                    break
                if filterReg.chromStart < snp.chromStart:
                    j += 1
                i +=1
                if i == len(filter):
                    break
                filterReg = filter[i]

        if inrange and snp.type == 'mnp' and not infilter:
            #sys.stderr.write("mnp\n")
            subsnps = splitDbSnp(snp)
            nummnp += len(subsnps)
            for s in subsnps:
                if isInRange(s.chromStart, start, end):
                    if len(snps) == 0:
                        snps.append( s )
                    else:
                        prevsnp = snps[ len(snps) - 1 ]
                        if prevsnp.chromStart == s.chromStart and prevsnp.chromEnd == s.chromEnd and prevsnp.chrom == s.chrom and prevsnp.observed == s.observed:
                            continue
                        snps.append(s)
        elif snp.isSnp and inrange and not infilter:
            if len(snps) == 0:
                snps.append( snp )
            else:
                prevsnp = snps[ len(snps) - 1 ]
                if prevsnp.chromStart == snp.chromStart and prevsnp.chromEnd == snp.chromEnd and prevsnp.chrom == snp.chrom and prevsnp.observed == snp.observed:
                    continue
                snps.append(snp)
    sys.stdout.write("TotalSnps\t%d\n" % (len(snps)))
    sys.stderr.write("numMNP: %d\n" %(nummnp))

    return snps

def readPgSnp(file, start, end, filter):
    f = open(file, 'r')
    sample2snps = {}
    i = 0
    j = 0
    
    for line in f:
        if re.search(line, 'chromStart'):
            continue
        snp = PgSnp(line)

        i = j
        inFilter = False
        if i < len(filter):
            filterReg = filter[i]
            while filterReg.chromStart <= snp.chromStart:
                if filterReg.chromStart <= snp.chromStart and snp.chromStart < filterReg.chromEnd:
                    inFilter = True
                    break
                if filterReg.chromStart < snp.chromStart:
                    j += 1
                i +=1
                if i == len(filter):
                    break
                filterReg = filter[i]

        #if not snp.isSnp or not isInRange(snp.chromStart, start, end) or isInFilter(snp.chromStart, filter):
        if not snp.isSnp or not isInRange(snp.chromStart, start, end) or inFilter:
            continue
        if snp.sample not in sample2snps:
            sample2snps[snp.sample] = [snp]
        else:
            sample2snps[snp.sample].append(snp)
    f.close()
    return sample2snps

def readRefSnps(file, filteredSamples, start, end, filter, longnames):
    xmltree = ET.parse( file )
    root = xmltree.getroot()
    snps = {}
    samples = []
    for sample in root.findall( 'statsForSample' ):
        name = sample.attrib['sampleName']
        totalCalls = int(sample.attrib['totalCalls'])
        if name == 'ROOT' or name == '' or name in filteredSamples:
            continue
        samples.append( (name, totalCalls) )
        snps[name] = []
        ref = sample.attrib['referenceName']
        #totalSnps = int( sample.attrib['totalErrors'] )
        sites = []
        if sample.text != None:
            sites = sample.text.strip().split('\n')
        for site in sites:#each snp detected by cactus ref, check to see if there is one in dbsnp
            if site != '':
               snp = SnpSite(site, name, ref, longnames)
               if isInRange(snp.refstart, start, end) and not isInFilter(snp.refstart, filter):
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
   
def calcOverlapSnps(snps, reportedSnps, isPgSnp, fpFh):#reportedSnps can be dbSnps or pgSnps
    refTotal = len(snps)
    totalSnps = len(reportedSnps)
    tpPos = 0 #pass position check (coordinate check) 
    tp = 0 #pass both coordinate and allele check
    currindex = 0

    for s in snps:
        flag = False
        flagTp = False
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
                    flagTp = True
                    break
                #else:
                #    sys.stderr.write("%s\t%d\t%d\tourref: %s\tourAllele: %s\tncbi: %s\tucsc: %s\tobserved: %s\n" %(s.sampleName, s.refstart, dbs, s.refallele, s.allele, reportedSnps[i].refNCBI, reportedSnps[i].refUCSC, ','.join(reportedSnps[i].observed) ))
            else: #s.refstart < dbs
                #print "\tGone too far, stop: %d, %d" %(s.refstart, dbs)
                break
        if flag == True: #pass position check
            tpPos += 1
        if not flagTp and fpFh != None:
            #fpFh.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n" %(s.sampleName, s.refchrom, s.refstart, s.refstart+1, s.allele, 1, 0, 0))
            fpFh.write("%s\t%s\t%d\t%d\t%s\t%s\t\t%s\t%d\n" %(s.sampleName, s.refchrom, s.refstart, s.refstart+1, s.allele.upper(), s.refallele.upper(), s.name, s.start))

    return tp, tpPos

def hackFalseNeg( sample, numSnpsConfirmed ):
    sample2totalSnps = {"cox": 15967, "qbl": 15282, "ssto": 14982, "apd": 4230, 'dbb': 14255, 'mann': 12102, 'mcf': 10790}
    if sample in sample2totalSnps:
        total = sample2totalSnps[sample]
        if total < numSnpsConfirmed:
            return 0, 0
        else:
            fn =  total - numSnpsConfirmed
            return fn, fn*100.0/total
    return -1, -1

def pgSnpOverlap(sample, snps, sample2snps, refTotal, tp): 
    #check against snps that were reported specifically for the sample
    if sample in sample2snps:
        pgsnps = sample2snps[sample]
        pgsnps.sort()
        pgTp, pgTpPos = calcOverlapSnps(snps, pgsnps, True, None)
        pgTotal = len(pgsnps)
        fn = pgTotal - pgTp
        if refTotal > 0 and pgTotal > 0:
            #sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
            sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 100.0*pgTpPos/refTotal, pgTp, 100.0*pgTp/refTotal, fn, 100.0*fn/pgTotal))
        else:
            sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 0.00, pgTp, 0.00, fn, 0.00))
    else:
        fn, fnPercentage = hackFalseNeg( sample, tp )
        if fn >= 0:
            sys.stdout.write("\t\t\t\t\t\t%d\t%.2f" %(fn, fnPercentage))

def getStats(dbsnps, refsnps, samples, sample2snps, sample2snpsPileup, falsePosFile):
    dbsnps.sort()

    fpFh = open(falsePosFile, "w")
    totalSnps = len(dbsnps)
    sys.stdout.write("Sample\tTotalCalledSnps\ttpPos\tPercentageTpPos\tTP\tPercentageTP\tFP\tPercentageFP\tsampleSnps\tsampleTpPos\tPercentageSampleTpPos\tsampleTP\tPercentageSampleTP\tsampleFN\tPercentageSampleFN\t")
    sys.stdout.write("pileupSnps\tpileupTpPos\tPercentagePileupTpPos\tpileupTP\tPercentagePileupTP\tpileupFN\tPercentagePileupFN\tTotalBases\n")
    for (sample, totalbases) in samples:
        #sys.stderr.write("Sample %s, sorting...\n" %(sample))
        snps = sorted( refsnps[sample] )
        #sys.stderr.write("Done sorting, len %d\n" %(len(snps)))
        
        #Check against dbSnps
        tp, tpPos = calcOverlapSnps(snps, dbsnps, False, fpFh)
        refTotal = len(snps)
        fp = refTotal - tp 
        if refTotal > 0:
            #sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
            sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tpPos, 100.0*tpPos/refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
        else:
            sys.stdout.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(sample, refTotal, tpPos, 0.00, tp, 0.00, fp, 0.00))
        
        pgSnpOverlap(sample, snps, sample2snps, refTotal, tp)
        pgSnpOverlap(sample, snps, sample2snpsPileup, refTotal, tp)

        sys.stdout.write("\t%d\n" %(totalbases))
    fpFh.close()

def initOptions( parser ):
    #parser.add_option('-o', '--outfile', dest='outfile', default='', help='Output file. Default is stdout' )
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    parser.add_option('--pgSnp', dest='pgSnp', help="pgSnp file (snps found in each sample)")
    parser.add_option('--pileupSnp', dest='puSnp', help="snp calls using samtools mpileup to reads alignment file (bam), in same format with pgSnp file")
    parser.add_option('-s', '--startCoord', dest='startCoord', type = 'int', help='Snps upstream of this Start coordinate (base 0) will be ignored. If not specified, it is assumed that there is no upstream limit.')
    parser.add_option('-e', '--endCoord', dest='endCoord', type='int', help='Snps downstream of this End coordinate () will be ignored. If not specified, it is assumed that there is no downstream limit.')
    parser.add_option('-f', '--filter', dest='filter', help='File contain regions to ignore in the stats (format:chr\\tchromStart\\tchromEnd). Will ignore all snps lie within this region. Default=no filtering')
    parser.add_option('--falsePos', dest='fp', help='File to write FalsePositive Calls. Default = "snpsFP.txt"')
    parser.add_option('--seqLongNames', dest='seqNameFile', help='File that mapped the short sequence names to longer names. Default=%default. (This option is added because the browser display requires that the sequence names are alphanumeric while the longer names with sample.chr.start.len.... provides mapping info)')

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
    if options.puSnp and not os.path.exists(options.puSnp):
        parser.error("pileupSnp file %s does not exists." %(options.puSnp))
    if options.filter != None:
        if not os.path.exists( options.filter ):
            parser.error("Repeat file %s does not exist\n" %options.filter)
        options.filter = readFilter( options.filter )
    else:
        options.filter = []
    if options.fp == None:
        options.fp = "snpsFP.txt"
    options.seqLongNames = {}
    if options.seqNameFile != None:
        options.seqLongNames = readSeqNameFile(options.seqNameFile)

def main():
    usage = ('Usage: %prog [options] snpStats_*.xml snp134dump.txt')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    
    dbsnps = readDbSnps( args[1], options.startCoord, options.endCoord, options.filter )
    sample2snps = {}
    sample2snpsPileup = {}
    if options.pgSnp:
        sample2snps = readPgSnp(options.pgSnp, options.startCoord, options.endCoord, options.filter)
    if options.puSnp:
        sample2snpsPileup = readPgSnp(options.puSnp, options.startCoord, options.endCoord, options.filter)
    refsnps,samples = readRefSnps( args[0], options.filteredSamples, options.startCoord, options.endCoord, options.filter, options.seqLongNames )
    getStats( dbsnps, refsnps, samples, sample2snps, sample2snpsPileup, options.fp )

    #Delete dbfile, refdbfile ...

if __name__ == '__main__':
    main()

