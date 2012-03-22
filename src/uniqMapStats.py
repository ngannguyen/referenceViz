#!/usr/bin/env python

"""
nknguyen at soe ucsc edu
Nov 23 2011

indir/
    experiment/
        catusRef/
            illumina/
                mergeSorted.bam
                snpcount.txt
                paired/
                    snpcount.txt
                single/
                    snpcount.txt
            summaryStats.txt
            cactusRef.fa
    otherRefs/
        sample/
            hg19/
                same structure with CactusRef

Outdir/
    experiment (required-...)/
        hg19/
            sample(e.g:NA12878)/
                
        apd/
        cox/
        ...
"""

import os, sys, re, copy
from optparse import OptionParser
from math import *

#from numpy import *
#import libPlotting as libplot
#import matplotlib.pyplot as pyplot
#from matplotlib.font_manager import FontProperties

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import system
from sonLib.bioio import logger
from sonLib.bioio import setLogLevel
from sonLib.bioio import getTempDirectory

class Setup(Target):
    def __init__(self, options):
        Target.__init__(self, time=0.0025)
        self.options = options

    def run(self):
        exps = getExps(self.options.indir)
        refIndir = os.path.join(self.options.indir, "otherRefs")
        for exp in exps:
            samples = exps[exp]
            for sample in samples:
                refs = os.listdir( os.path.join(refIndir, sample) )
                for ref in refs:
                    outdir = os.path.join(self.options.outdir, exp, sample, ref)
                    system("mkdir -p %s" %outdir)
                    self.addChildTarget( RunExp(exp, ref, sample, outdir, self.options) )

class RunExp(Target):
    def __init__(self, exp, ref, sample, outdir, options):
        Target.__init__(self)
        self.exp = exp
        self.ref = ref
        self.sample = sample
        self.outdir = outdir
        self.options = options

    def run(self):
        expdir = os.path.join(self.options.indir, "_".join([self.exp, self.sample]))
        cactusRefDir = os.path.join( expdir, "cactusRef" )
        cactusRef = os.path.join(cactusRefDir, "cactusRef.fa")
        
        otherRefDir = os.path.join(self.options.indir, "otherRefs", self.sample, self.ref)
        otherRef = os.path.join(otherRefDir, self.ref)
        
        #Copy files to global dir:
        localTempDir = self.getLocalTempDir()
        crefLocalDir = os.path.join(localTempDir, "cref")
        system("mkdir -p %s" %crefLocalDir)
        system("scp -C %s %s" %(cactusRef, crefLocalDir))
        system("samtools faidx %s" % os.path.join(crefLocalDir, 'cactusRef.fa'))
        system("scp -C %s %s" %(os.path.join(cactusRefDir, "illumina", "mergeSorted.bam"), os.path.join(crefLocalDir, "CRef.bam")) )
        system("scp -C %s %s" %(os.path.join(cactusRefDir, "illumina", "merge.vcf"), crefLocalDir) )
        system("scp %s %s" %(os.path.join(cactusRefDir, "illumina", "snpcount.txt"), crefLocalDir) )

        orefLocalDir = os.path.join(localTempDir, 'oref')
        system("mkdir -p %s" %(orefLocalDir))
        system("scp -C %s %s" %(otherRef, orefLocalDir))
        system("samtools faidx %s" % os.path.join(orefLocalDir, '%s'%self.ref))
        system("scp -C %s %s" %(os.path.join(otherRefDir, "illumina", "mergeSorted.bam"), os.path.join(orefLocalDir, "%s.bam" %self.ref)) )
        system("scp -C %s %s" %(os.path.join(otherRefDir, "illumina", "merge.vcf"), orefLocalDir) )
        system("scp %s %s" %(os.path.join(otherRefDir, "illumina", "snpcount.txt"), orefLocalDir) )

        #Compare uniq versus multi mapping reads betwe`en cref and oref:
        increfU = os.path.join(self.outdir,  "%sM-CRefU_U.bam" %self.ref)
        inorefU = os.path.join(self.outdir,  "CRefM-%sU_U.bam" %self.ref)
        
        crefU = os.path.join(localTempDir, "%sM-CRefU_U.bam" %self.ref)
        orefU = os.path.join(localTempDir, "CRefM-%sU_U.bam" %self.ref)
        if os.path.exists(increfU) and os.path.exists(inorefU):
            system("scp -C %s %s" %(increfU, localTempDir))
            system("scp -C %s %s" %(inorefU, localTempDir))
        else:
            system("bam_cmpMapping.py %s %s %s" %( os.path.join(crefLocalDir, "CRef.bam"), os.path.join(orefLocalDir, "%s.bam" %self.ref), localTempDir))
            system("scp -C %s %s" %(crefU, self.outdir))
            system("scp -C %s %s" %(orefU, self.outdir))

        #====== Get number of reads that are uniquely mapped to one ref but not to the other:
        crefUcount = os.path.join(localTempDir, "crefUcount.txt") 
        system("samtools view %s | wc -l > %s" %(crefU, crefUcount))
        
        orefUcount = os.path.join(localTempDir, "orefUcount.txt") 
        system("samtools view %s | wc -l > %s" %(orefU, orefUcount))
       
        #=========== PILEUP ===========
        #Pileup for CRef:
        crefUsorted = os.path.join(localTempDir, "%sM-CRefU_U-sorted" %self.ref)
        system("samtools sort %s %s" %(crefU, crefUsorted))
        crefUpileup, crefUvcf, crefPileup = snpCalling( crefLocalDir, "%s.bam" % crefUsorted, os.path.join(crefLocalDir, "CRef.bam"), "cactusRef.fa")

        #Pileup for otherRef
        orefUsorted = os.path.join(localTempDir, "CRefM-%sU_U-sorted" %self.ref)
        system("samtools sort %s %s" %(orefU, orefUsorted))
        orefUpileup, orefUvcf, orefPileup = snpCalling( orefLocalDir, "%s.bam" % orefUsorted, os.path.join(orefLocalDir, "%s.bam" %self.ref), self.ref)

        #=========== SNP COUNTS ===================
        crefSnpcount, crefPileupCount, crefUsnpcount, crefUpileupCount = snpCount(crefLocalDir, crefUpileup, crefUvcf, crefPileup)
        orefSnpcount, orefPileupCount, orefUsnpcount, orefUpileupCount = snpCount(orefLocalDir, orefUpileup, orefUvcf, orefPileup)
        
        #=========
        crefSites = readPileup(crefUpileup)
        crefIntervals = covFilter(crefSites, 0)
        crefLendist = lenDist(crefIntervals, crefLocalDir)
        crefRepeat = repeat(crefIntervals, crefLocalDir)
        crefVcf = os.path.join(crefLocalDir, 'merge.vcf') 
        crefSnpcount2 = getSnpSubRegions(crefSites, crefVcf, crefLocalDir)

        orefSites = readPileup(orefUpileup)
        orefIntervals = covFilter(orefSites, 0)
        orefLendist = lenDist(orefIntervals, orefLocalDir)
        orefRepeat = repeat(orefIntervals, orefLocalDir)
        orefVcf = os.path.join(orefLocalDir, 'merge.vcf') 
        orefSnpcount2 = getSnpSubRegions(orefSites, orefVcf, orefLocalDir)
        
        #======== PRINT TO OUTPUT FILE =============
        crefout = printCount(crefUcount, crefPileupCount, crefSnpcount, crefUpileupCount, crefUsnpcount, crefSnpcount2, crefLocalDir)
        orefout = printCount(orefUcount, orefPileupCount, orefSnpcount, orefUpileupCount, orefUsnpcount, orefSnpcount2, orefLocalDir)

        #======== Copy outfiles to output directory =========
        crefOutdir = os.path.join(self.outdir, 'cref')
        system("mkdir -p %s" %crefOutdir)
        orefOutdir = os.path.join(self.outdir, self.ref)
        system("mkdir -p %s" %orefOutdir)

        system("cp %s %s" %(crefout, crefOutdir))
        system("cp %s %s" %(crefLendist, crefOutdir))
        system("cp %s %s" %(crefRepeat, crefOutdir))
        system("cp %s %s" %(orefout, orefOutdir))
        system("cp %s %s" %(orefLendist, orefOutdir))
        system("cp %s %s" %(orefRepeat, orefOutdir))

        #======== Test with dbSNP ==========
        dbsnps = readDbSNPs(self.options.dbsnp, self.options.startCoord, self.options.endCoord) 
        dbsnpFile = getDbSnpStats(dbsnps, orefSites, self.options.startCoord, self.options.endCoord, orefLocalDir)
        system("cp %s %s" %(dbsnpFile, orefOutdir))

        #self.addChildTarget( UniqHighCov(cactusRef, cactusRefDir, os.path.join(self.outdir, 'cref'), crefCovCutoff) )
        #self.addChildTarget( UniqHighCov(otherRef, otherRefDir, os.path.join(self.outdir, self.ref), orefCovCutoff) )
        
        #self.setFollowOnTarget()

class UniqHighCov(Target):
    def __init__(self, ref, indir, outdir, covCutoff):
        Target.__init__(self)
        self.ref = ref
        self.indir = indir
        self.outdir = outdir
        self.covCutoff = covCutoff

    def run(self):
        system("mkdir -p %s" %(self.outdir))
        localTempDir = self.getLocalTempDir()
        bam = os.path.join(localTempDir, "mergeSorted.bam")
        localRef = os.path.join(localTempDir, os.path.basename(self.ref) )
        system( "scp -C %s %s" %(os.path.join(self.indir, "illumina", "mergeSorted.bam"), bam ))
        system( "scp -C %s %s" %(self.ref, localTempDir) )
        system( "samtools faidx %s" %( localRef ) )

        #Get Uniquely mapped reads:
        uniqbam = os.path.join(localTempDir, "uniq.bam")
        uniqSorted = os.path.join(localTempDir, "uniqSorted")
        if os.path.exists( os.path.join(self.outdir, "uniqSorted.bam") ):
            system("scp -C %s %s" %(os.path.join(self.outdir, "uniqSorted.bam"), localTempDir))
        else:
            system( "getBamUniqMappedReads.py < %s > %s" %(bam, uniqbam) )
            #Sort
            system( "samtools sort %s %s" %(uniqbam, uniqSorted) )
            system( "scp -C %s.bam %s" %(uniqSorted, self.outdir))

        #Pileup:
        pileup = os.path.join(localTempDir, "pileup.txt")
        system("samtools mpileup %s.bam -f %s > %s" %(uniqSorted, localRef, pileup))
        system("scp -C %s %s" %(pileup, self.outdir))
        
        #Analyze pileup:
        sites = readPileup(pileup)
        intervals = covFilter(sites, self.covCutoff)
        lenDist(intervals, self.outdir)
        repeat(intervals, self.outdir)
        
class Site():
    def __init__(self, line):
        items = line.strip().split()
        if len(items) < 4:
            sys.stderr.write("Invalid pileup entry:\n%s\n" %line)
            sys.exit(1)
        self.ref = items[0]
        self.pos = int(items[1])
        self.base = items[2]
        self.cov = int(items[3])

def getSnpSubRegions(sites, vcf, dir):
    #sites = readPileup(upileup)
    siteList = os.path.join(dir, "uniqSites.txt")
    f = open(siteList, 'w')
    for site in sites:
        f.write("%s\t%d\n" %(site.ref, site.pos))
    f.close()
    outfile = os.path.join(dir, "uniqSnpcount2.txt")
    system("bcftools view -S -l %s %s | grep -vc '^#' > %s" %(siteList, vcf, outfile))
    #count = getCount(outfile) 
    return outfile
        
def lenDist(intervals, outdir):
    #intervals = covFilter(sites, self.covCutoff)
    len2count = getLenDist(intervals)
    totalCount = sum( [ len2count[l] for l in len2count ] )
    lenDistOut = os.path.join(outdir, "lenDist.txt")
    if totalCount > 0:
        f1 = open(lenDistOut, 'w')
        for l in sorted( len2count.keys() ):
            f1.write("%d\t%d\t%.3f\n" %(l, len2count[l], 100.0*len2count[l]/totalCount ))
        f1.close()
    return lenDistOut

def repeat(intervals, outdir):
    repeat, total = getRepeat(intervals)
    rpOut = os.path.join(outdir, "repeat.txt")
    if total > 0:
        f2 = open(rpOut, 'w')
        f2.write("%d\t%d\t%.3f\n" %(repeat, total, 100.0*repeat/total))
        f2.close()
    return rpOut

def snpCalling(dir, uniqBam, bam, ref):
    reffile = os.path.join(dir, ref)
    ##pileup of reads the are uniquely mapped to Cref, but not to hg19
    upileup = os.path.join(dir, "pileupU.txt")
    system("samtools mpileup %s -f %s > %s" %(uniqBam, reffile, upileup))
    ##SNP calls for reads the are uniquely mapped to ref, but not to otherRef
    uvcf = os.path.join(dir, "unique.vcf")
    system("samtools mpileup %s -f %s -g | bcftools view -v -I - > %s" %(uniqBam, reffile, uvcf))
    ##pileup of all reads
    pileup = os.path.join(dir, "mergePileup.txt")
    system("samtools mpileup %s -f %s > %s" %( bam, reffile, pileup))
    #crefCovCutoff = getCovCutoff(crefPileup, self.options.covCutoff)
    #Get SNP calls for regions that uniquely-mapped reads mapped to:
    #crefSnp = getSnpSubRegions(crefPileup, os.path.join(crefLocalDir, "CRef.bam"), os.path.join(crefLocalDir, "cactusRef.fa"))
    return upileup, uvcf, pileup
    
def snpCount(dir, upileup, uvcf, pileup):
    snpcount = os.path.join(dir, "snpcount.txt")
    
    pileupcount = os.path.join(dir, "mergePileupCount.txt")
    system("grep -vc '^#' %s > %s" %(pileup, pileupcount))

    usnpcount = os.path.join(dir, "uniqSnpcount.txt")
    system("grep -vc '^#' %s > %s" %(uvcf, usnpcount))

    upileupcount = os.path.join(dir, "uniqPileupCount.txt")
    system("grep -vc '^#' %s > %s" %(upileup, upileupcount))
    return snpcount, pileupcount, usnpcount, upileupcount

def getCovCutoff( pileupFile, percentile ):
    f = open(pileupFile, 'r')
    covdist = {} #key = cov, val = count
    for line in f:
        if line[0] == '#':
            continue
        items = line.strip().split()
        if len(items) < 4:
            continue
        cov = int(items[3])
        if cov in covdist:
            covdist[cov] +=1
        else:
            covdist[cov] = 1
    total = sum([covdist[c] for c in covdist])
    if total == 0:
        return 0
    percentileSum = float(total)*percentile
    count = 0
    prevcov = 0
    for cov in sorted( covdist.keys() ):
        if count + covdist[cov] > percentileSum:
            return prevcov
        else:
            count += covdist[cov]
            prevcov = cov
    return prevcov

def covFilter(sites, cutoff):
    intervals = [] #list of interval of sites with cov >= cutoff
    currlist = []
    for site in sites:
        if site.cov < cutoff:
            if len(currlist) > 0:
                intervals.append(currlist)
                currlist = []
            continue
        else:
            if len(currlist) == 0:
                currlist.append(site)
            else:
                prevSite = currlist[ len(currlist) -1 ]
                if prevSite.pos ==  site.pos -1 :
                    currlist.append(site)
                else:
                    intervals.append(currlist)
                    currlist = [site]
    if len(currlist) > 0:
        intervals.append(currlist)
    return intervals

def getRepeat(intervals):
    repeat = 0
    total = 0
    for interval in intervals:
        for site in interval:
            total +=1
            if site.base == site.base.lower():
                repeat += 1
    return repeat, total

def getLenDist(intervals):
    len2count = {} #key = len, val = count
    for interval in intervals:
        start = interval[0].pos
        end = interval[ len(interval) -1 ].pos 
        length = end - start + 1
        if length in len2count:
            len2count[length] +=1
        else:
            len2count[length] = 1
    return len2count

def readPileup(file):
    f = open(file, 'r')
    sites = []
    for line in f:
        if line[0] != '#':
            sites.append(Site(line))
    return sites

#crefout = printCount(crefPileupCount, crefSnpcount, crefUpileupCount, crefUsnpcount, crefSnpcount2, crefLocalDir)
def printCount(uReads, pileup, snp, uniqPileup, uniqSnp, uniqSnp2, outdir):
    outfile = os.path.join(outdir, 'uniqMapStats.txt')
    f = open(outfile, 'w')
    uReadsCount = getCount(uReads)
    pileupCount = getCount(pileup)
    snpCount = getCount(snp)
    uniqPileupCount = getCount(uniqPileup)
    uniqSnpCount = getCount(uniqSnp)
    uniqSnpCount2 = getCount(uniqSnp2)
    f.write("#Number of reads that are uniquely mapped to this ref, but not the other ref: %d\n" %uReadsCount)
    f.write("#Category\tSnpCount\tTotalBases\tSnpPerBase\n")
    if pileupCount > 0:
        f.write("All reads\t%d\t%d\t%.3f\n" %(snpCount, pileupCount, 100.0*snpCount/pileupCount))
    else:
        f.write("All reads\t%d\t%d\t%.3f\n" %(0, 0, 0.0))
    if uniqPileupCount > 0:
        f.write("All reads-uniqM Reg\t%d\t%d\t%.3f\n" %(uniqSnpCount2, uniqPileupCount, 100.0*uniqSnpCount2/uniqPileupCount))
        f.write("UniqM reads\t%d\t%d\t%.3f\n" %(uniqSnpCount, uniqPileupCount, 100.0*uniqSnpCount/uniqPileupCount))
    else:
        f.write("All reads-uniqM reg\t%d\t%d\t%.3f\n" %(0, 0, 0.0))
        f.write("UniqM reads\t%d\t%d\t%.3f\n" %(0, 0, 0.0))
    return outfile

def getCount(file):
    f = open(file, 'r')
    line = f.readline()
    return int(line.strip())

def getExps(indir):
    dirs = os.listdir(indir)
    exps = {} #key = exp, val = [samples]
    for d in dirs:
        if re.search('required', d):
            items = d.split('_')
            if len(items) == 3:
                sample = items[2]
                exp = "_".join(items[:2])
                if exp not in exps:
                    exps[exp] = [sample]
                else:
                    exps[exp].append(sample)
    return exps

#==================== dbSNP ===============================
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
            self.observed = [ reverse(a) for a in self.observed] 
        #Check if indel is snp:
        self.isSnp = False
        for o in self.observed:
            if len(o) == self.chromEnd - self.chromStart and o!= '-' and o != self.refUCSC:
                self.isSnp = True
                break
        
        self.molType = items[9]
        self.type = items[10]
        self.func = items[14]
        self.locType = items[15]
        self.exceptions = items[17]
        self.alleles = items[21].lower().rstrip(',').split(',')
        self.alleleFreqs = {}
        if items[23] != "":
            alleleFreqs = [ float(freq) for freq in items[23].rstrip(',').split(',') ]
            for i, allele in enumerate(self.alleles):
                self.alleleFreqs[ allele ] = alleleFreqs[i]
        if self.strand == '-':
            d = {'a':'t', 't':'a', 'c':'g', 'g':'c'}
            newalleles = []
            for a in self.alleles:
                if a in d:
                    newalleles.append( d[a] )
                else:
                    newalleles.append( a )
            self.alleles = newalleles

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

def splitDbsnp(snp):
    snps = []
    if len(snp.observed) == 0:
        return snps
    snpLen = len(snp.observed[0])
    for allele in snp.observed:
        if len(allele) != snpLen:
            return snps
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

def readDbSNPs(file, start, end):
    f = open(file, 'r')
    snps = []

    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        inrange = isInRange(snp.chromStart, start, end)
        if inrange and snp.type == 'mnp':
            subsnps = splitDbsnp(snp)
            for s in subsnps:
                if len(snps) == 0:
                    snps.append(s)
                else:
                    prevsnp = snps[ len(snps) - 1 ]
                    #if prevsnp.chromStart == s.chromStart and prevsnp.chromEnd == s.chromEnd and prevsnp.chrom == s.chrom and prevsnp.observed == s.observed:
                    if prevsnp.chromStart == s.chromStart and prevsnp.chromEnd == s.chromEnd and prevsnp.chrom == s.chrom:
                        continue
                    snps.append(s)
        elif inrange and snp.isSnp:
            if len(snps) == 0:
                snps.append(snp)
            else:
                prevsnp = snps[ len(snps) - 1 ]
                if prevsnp.chromStart == snp.chromStart and prevsnp.chromEnd == snp.chromEnd and prevsnp.chrom == snp.chrom:
                    continue
                snps.append(snp)

    return snps

def getDbSnpStats(dbsnps, sites, start, end, outdir):
    totalbases = end - start
    totalSnps = len(dbsnps)
    if totalbases <=0 :
        exit(1)
    overallrate = float(totalSnps)/totalbases
    #uniquely mapped regions:
    uniqbases = len(sites)
    i = 0
    uniqsnps = 0
    for snp in dbsnps:
        while i < len(sites) and sites[i].pos + start <= snp.chromStart:
            if sites[i].pos + start == snp.chromStart:
                uniqsnps += 1
                break
            i += 1
    uniqrate = float(uniqsnps)/uniqbases
    ratio = 0.0
    if overallrate > 0:
        ratio = uniqrate/overallrate
    
    #Print to output file
    outfile = os.path.join(outdir, 'uniqMap-dbsnp.txt')
    f = open(outfile, 'w')
    f.write("Category\tTotalbases\tTotalSNPs\tRate\n")
    f.write("Overall\t%d\t%d\t%f\n" %(totalbases, totalSnps, overallrate))
    f.write("UniquelyMapped\t%d\t%d\t%f\n" %(uniqbases, uniqsnps, uniqrate))
    f.write("Ratio\t%f\n" %(ratio))

    return outfile

#=================== END dbSNP ============================

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory. Required argument.')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory. Default = ./uniqMapping')
    parser.add_option('-c', '--covCutoff', dest='covCutoff', type='int', default = 0.9, help='Coverage cutoff. Default = 0.9' )
    parser.add_option('-d', '--dbSNP', dest='dbsnp', help='dbSNPfile', default='/hive/users/nknguyen/reconGit/referenceScripts/data/snp134.txt')
    parser.add_option('-s', '--startCoord', dest='startCoord', type = 'int', default=28477754, help='Snps upstream of this Start coordinate (base 0) will be ignored. If not specified, it is assumed that there is no upstream limit.')
    parser.add_option('-e', '--endCoord', dest='endCoord', type='int', default=33448354, help='Snps upstream of this Start coordinate (base 0) will be ignored. If not specified, it is assumed that there is no downstream limit.')
 
    #parser.add_option()

def checkOptions( args, options, parser ):
    if options.indir == None:
        parser.error('No input directory was given.\n')
    if not os.path.exists( options.indir ):
        parser.error('Input directory does not exist.\n')
    if options.outdir == None:
        options.outdir = os.path.join( os.getcwd(), "mapStats" )
    #options.refs = options.refs.split(',')
    #system("mkdir -p %s" % options.outdir)

def main():
    usage = "uniqMapStats.py [options]"
    parser = OptionParser(usage = usage)
    Stack.addJobTreeOptions(parser)

    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs" %i)

if __name__ == '__main__':
    from referenceViz.src.uniqMapStats import *
    main()
