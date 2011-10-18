#!/usr/bin/env python

"""
Compare Snps between each sample and the reference with the set of Snps in dbSnps
"""

import os, sys, re
#import sqlite3
from optparse import OptionParser
import xml.etree.ElementTree as ET
from sonLib.bioio import system

class Filter():
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 3:
            sys.stderr.write("Filter region has wrong format: %s\n" %line)
        self.chrom = items[0]
        self.chromStart = int(items[1])
        self.chromEnd = int(items[2])

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
        self.alleleCount = int(items[4])
        if self.alleleCount == 0:
            self.alleleSize = 0
        else:
            self.alleleSize = len(self.alleles[0])
        self.isIndel = ''
        if re.search('-', self.name):
            if self.chromEnd == self.chromStart: #insertion
                self.isIndel = 'insertion'
            else:
                self.isIndel = 'deletion'


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
        if len(items) < 16:
            sys.stderr.write('cactus indel record does not have enough fields\n')
            sys.exit(1)
        self.name = items[9]
        self.start = int( items[11] )
        self.length = int(items[13])
        self.strand = items[15]
        #self.allele = items[2].lower()
        self.ref = items[1]
        self.refstart = int( items[3] ) 
        self.reflength = int(items[5])
        self.refstrand = items[7]
        if self.reflength == 0 and self.length == 0:
            sys.stderr.write('Not an insertion or a deletion: length1 == length2 ==0: %s\t%d\n' %(self.sampleName, self.refstart))

        refitems = self.ref.split('.')
        if len(refitems) <=2:
            sys.stderr.write('Reference sequence must provide chromosome info: i.e: species.chr\n')
            sys.exit(1)
        self.refchrom = refitems[1]
        
        #TEMP:
        if self.refstrand == 0:
            sys.stderr.write("Cactus indels are on negative strand\n")
            sys.exit(1)

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
        #self.refallele = items[5].lower()
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

def checkSizeCutoff(snp, cutoff):
    #Note, only check for 'deletion' and insertion' cases. In case of in-del, don't check
    if snp.type == 'deletion' and snp.chromEnd - snp.chromStart > cutoff:
        return False
    elif snp.type == 'insertion':
        flag = False
        for allele in snp.observed:
            if allele != '-' and len(allele) <= cutoff:
                flag = True
                break
        return flag
    return True

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

def readDbSnps(file, cutoff, start, end, filter):
    f = open(file, 'r')
    #indels = []
    insertions = []
    deletions = []
    numins = 0
    numdels = 0
    numindels = 0

    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        inrange = isInRange(snp.chromStart, start, end)
        infilter = isInFilter( snp.chromStart, filter )
        if not inrange or infilter:
            continue
        if snp.type == 'in-del':
            #indels.append(snp)
            numindels +=1
            insertions.append(snp)
            deletions.append(snp)
        elif snp.type == 'insertion' and checkSizeCutoff(snp, cutoff):
            numins +=1
            insertions.append(snp)
        elif snp.type == 'deletion' and checkSizeCutoff(snp, cutoff):
            numdels += 1
            deletions.append(snp)
    
    sys.stdout.write("Insertions\t%d\tDeletions\t%d\tIn-dels\t%d\n" %(numins, numdels, numindels))

    #return indels, insertions, deletions
    return insertions, deletions

def readPgSnp(file, cutoff, start, end, filter):
    f = open(file, 'r')
    sample2ins = {}
    sample2dels = {}
    for line in f:
        if re.search(line, 'chromStart'):
            continue
        snp = PgSnp(line)
        if snp.isIndel != ""  and snp.alleleSize <= cutoff and isInRange(snp.chromStart, start, end) and not isInFilter(snp.chromStart, filter):
            if snp.isIndel == 'insertion':
                if snp.sample not in sample2ins:
                    sample2ins[snp.sample] = [snp]
                else:
                    sample2ins[snp.sample].append(snp)
            elif snp.isIndel == 'deletion':
                if snp.sample not in sample2dels:
                    sample2dels[snp.sample] = [snp]
                else:
                    sample2dels[snp.sample].append(snp)
                
    f.close()
    return sample2ins, sample2dels

def readRefSnps(file, filteredSamples, cutoff, start, end, filter):
    xmltree = ET.parse( file )
    root = xmltree.getroot()
    ins = {}
    dels = {}
    samples = []
    for sample in root.findall( 'statsForSample' ):
        name = sample.attrib['sampleName']
        if name == 'ROOT' or name == '' or name in filteredSamples:
            continue
        samples.append(name)
        ins[name] = []
        dels[name] = []
        ref = sample.attrib['referenceName']
        #totalSnps = int( sample.attrib['totalErrors'] )
        sites = sample.text.strip().split('\n')
        for site in sites:#each snp detected by cactus ref, check to see if there is one in dbsnp
            if site != '':
                snp = SnpSite(site, name, ref)
                if snp.length <= cutoff and snp.reflength <= cutoff and isInRange(snp.refstart, start, end) and not isInFilter(snp.refstart, filter):
                    if snp.reflength > 0 and snp.length == 0:
                        dels[name].append(snp)
                    elif snp.reflength == 0 and snp.length > 0:
                        ins[name].append(snp)
                    else:
                        ins[name].append(snp)
                        dels[name].append(snp)

    return ins, dels, samples

def checkAlleles(refsnp, dbsnp, type):
    if type == 'insertion':
        reflen = refsnp.length
    elif type == 'deletion':
        reflen = refsnp.reflength
    else:
        sys.stderr.write('Unknown type %s. Expected value is "deletion" or "insertion".\n')
        sys.exit(1)
    #reflen = abs(refsnp.length - refsnp.reflength) #our indel length
    dblen = dbsnp.chromEnd - dbsnp.chromStart
    if dbsnp.type == 'deletion':
        if reflen == dblen:
            return True
        else:
            return False
    elif dbsnp.type == 'insertion' or refsnp.length > 0:#insertion
        for allele in dbsnp.observed:
            if allele != '-' and len(allele) == reflen:
                return True
        return False
    elif refsnp.reflength > 0: #in-del, deletion
        return True
    return True

def checkPgAlleles(refsnp, pgsnp, type):
    if type == 'insertion':
        reflen = refsnp.length
    elif type == 'deletion':
        reflen = refsnp.reflength
    else:
        sys.stderr.write('Unknown type %s. Expected value is "deletion" or "insertion".\n')
        sys.exit(1)
    #reflen = abs(refsnp.length - refsnp.reflength)
    if reflen == pgsnp.alleleSize:
        return True
    return False

def calcDbIndelOverlap(snps, dbsnps, wobble, isPgSnp, type):#reportedSnps can be dbSnps or pgSnps
    refTotal = len(snps)
    totalSnps = len(dbsnps)
    tp = 0 
    tpPos = 0
    currindex = 0

    #dumpfh = open("dump.txt", "a")

    for s in snps: #each indel
        flag = False
        #tpflag = False
        #dbs = 0
        #dbe = 0
        #for i in xrange( totalSnps ):
        for i in xrange(currindex, totalSnps):
            dbs = dbsnps[i].chromStart
            dbe = dbsnps[i].chromEnd

            if dbs < s.refstart - wobble:
                currindex = i
            #elif flag == False:
            #    currindex = i
            if dbs < s.refstart - wobble:
                #print "\tnot there yet: %d < %d" %(s.refstart, dbs)
                continue
            elif  s.refstart - wobble <= dbs and dbs <= s.refstart + wobble: # < dbe
                flag = True
                #print "\tfound it! %d, %d" %(s.refstart, dbs)
                if (isPgSnp and checkPgAlleles(s, dbsnps[i], type) ) or (not isPgSnp and checkAlleles(s, dbsnps[i], type)):
                    tp += 1
                    #tpflag = True
                    #if s.sampleName != 'panTro3':
                    #    dumpfh.write("%d\t%s\trefStart: %d,%d\tStart: %d,%d\tdbIndels, start: %d, length: %d\n" %(s.refstart - dbs, s.sampleName, s.refstart, s.reflength, s.start, s.length, dbs, dbe - dbs))
                    break
                #else:
                #    sys.stderr.write("%s\t%d\t%d\tourref: %s\tourAllele: %s\tncbi: %s\tucsc: %s\tobserved: %s\n" %(s.sampleName, s.refstart, dbs, s.refallele, s.allele, reportedSnps[i].refNCBI, reportedSnps[i].refUCSC, ','.join(reportedSnps[i].observed) ))
            else: #s.refstart + wobble < dbs
                #print "\tGone too far, stop: %d, %d" %(s.refstart, dbs)
                break
        if flag == True: #pass position check
            tpPos += 1
        #else :
        #    dumpfh.write("%d\t%s\trefStart: %d,%d\tStart: %d,%d\tdbIndels, start: %d, length: %d\n" %(s.refstart - dbs, s.sampleName, s.refstart, s.reflength, s.start, s.length, dbs, dbe - dbs))

    return tp, tpPos

def getStats(dbsnps, refsnps, samples, sample2snps, wobble, type):
    dbsnps.sort()
    #dbindels.sort()

    totalSnps = len(dbsnps)
    for sample in samples:
        snps = sorted( refsnps[sample] )
        
        #Check against dbSnps (Insertions or deletions)
        tp, tpPos = calcDbIndelOverlap(snps, dbsnps, wobble, False, type)
        #tp1, tpPos1 = calcDbIndelOverlap(snps, dbindels, wobble, False)
        #tp += tp1
        #tpPos += tpPos1
        refTotal = len(snps)
        fp = refTotal - tp 
        if refTotal > 0:
            sys.stdout.write("%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(type, sample, refTotal, tpPos, 100.0*tpPos/refTotal, tp, 100.0*tp/refTotal, fp, 100.0*fp/refTotal))
        else:
            sys.stdout.write("%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(type, sample, refTotal, tpPos, 0.00, tp, 0.00, fp, 0.00))

        #check against snps that were reported specifically for the sample
        if sample in sample2snps:
            pgsnps = sample2snps[sample]
            pgsnps.sort()
            pgTp, pgTpPos = calcDbIndelOverlap(snps, pgsnps, wobble, True, type)
            pgTotal = len(pgsnps)
            fn = pgTotal - pgTp
            if refTotal > 0 and pgTotal > 0:
                sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 100.0*pgTpPos/refTotal, pgTp, 100.0*pgTp/refTotal, fn, 100.0*fn/pgTotal))
            else:
                sys.stdout.write("\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f" %(pgTotal, pgTpPos, 0.00, tp, 0.00, fn, 0.00))

        sys.stdout.write("\n")


def initOptions( parser ):
    #parser.add_option('-o', '--outfile', dest='outfile', default='', help='Output file. Default is stdout' )
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    parser.add_option('--pgSnp', dest='pgSnp', help="pgSnp file (snps found in each sample)")
    parser.add_option('-w', dest='wobble', type='int', default=0, help="Default = 0")
    parser.add_option('-c', '--cutoff', dest='cutoff', type='int', default = 10, help='Cutoff size of indels to be included in the analysis. Default = 10')
    parser.add_option('-s', '--startCoord', dest='startCoord', type = 'int', help='Snps upstream of this Start coordinate (base 0) will be ignored. If not specified, it is assumed that there is no upstream limit.')
    parser.add_option('-e', '--endCoord', dest='endCoord', type='int', help='Snps upstream of this Start coordinate (base 0) will be ignored. If not specified, it is assumed that there is no downstream limit.')
    parser.add_option('-f', '--filter', dest='filter', help='File contain regions to ignore in the stats (format:chr\\tchromStart\\tchromEnd). Will ignore all snps lie within this region. Default=no filtering')

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
    if options.filter != None:
        if not os.path.exists( options.filter ):
            parser.error("Repeat file %s does not exist\n" %options.filter)
        options.filter = readFilter( options.filter )
    else:
        options.filter = []

def main():
    usage = ('Usage: %prog [options] pathStats_*.xml snp134dump.txt')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    
    dbinsertions, dbdeletions = readDbSnps( args[1], options.cutoff, options.startCoord, options.endCoord, options.filter )
    sample2ins = {}
    sample2dels = {}
    if options.pgSnp:
        sample2ins, sample2dels = readPgSnp(options.pgSnp, options.cutoff, options.startCoord, options.endCoord, options.filter)
    refins, refdels,samples = readRefSnps( args[0], options.filteredSamples, options.cutoff, options.startCoord, options.endCoord, options.filter )
    sys.stdout.write("Type\tSample\tTotalCalled\ttpPos\tPercentageTpPos\tTP\tPercentageTP\tsampleSnps\tsampleTpPos\tPercentageSampleTpPos\tsampleTP\tPercentageSampleTP\tsampleFN\tPercentageSampleFN\n")   
    #sys.stdout.write("Insertions:\n")
    getStats( dbinsertions, refins, samples, sample2ins, options.wobble, 'insertion')
    #sys.stdout.write("\nDeletions:\n")
    getStats( dbdeletions, refdels, samples, sample2dels, options.wobble, 'deletion' )

    #Delete dbfile, refdbfile ...

if __name__ == '__main__':
    main()

