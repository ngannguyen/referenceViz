#!/usr/bin/env python

"""
nknguyen at soe ucsc edu
Sep 12 2011
Input: 
Output:
"""
import os, sys, re
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

class Setup(Target):
    """
    Setup mapping runs
    """
    def __init__(self, options):
        #Target.__init__(self, time=0.00025)
        Target.__init__(self)
        self.options = options

    def run(self):
        setLogLevel("DEBUG")
        options = self.options
        system("mkdir -p %s" %options.outdir)

        #Run bwa index:
        refpath = os.path.join( options.outdir, "bwaRef" )
        refpath = os.path.join( refpath, "bwaRef")
        if options.bwaRefIndex != None:
            refpath = options.bwaRefIndex
        else:
            system("mkdir -p %s" %refpath)
            system( "bwa index -p %s %s" %(refpath , options.ref) )
            logger.info("Bwa indexed reference %s at with prefix %s\n" %(options.ref, refpath) )
        self.addChildTarget( Paired(options, refpath, "illumina", self.options.iMaxInsertSize) )
        self.addChildTarget( Single(options, refpath, "illumina") )
        #GOt error running ssaha2 for 454paired (required 2 reads of the same pair to have the same length), so for now, use bwa
        self.addChildTarget( Paired(options, refpath, "ls454", self.options.ls454MaxInsertSize) )
        #bwa sw (default setting) sometimes return hits with low identity percentage, as well as return chimera hits (so one read can potentially have 2, 3 hits)
        self.addChildTarget( Single(options, refpath, "ls454") )

        #Run bwasw index (for ls454 single)
        #refpath = os.path.join( options.outdir, "bwaswRef" )
        #system("mkdir -p %s" %refpath)
        #refpath = os.path.join( refpath, "bwaswRef")
        #system( "bwa index -p %s -a bwtsw %s" %(refpath , options.ref) )
        #logger.info("Bwasw indexed reference %s at with prefix %s\n" %(options.ref, refpath) )
        #self.addChildTarget( Ls454Single(options, refpath) )
        
        #Run ssaha2 index:
        #refpath = os.path.join( options.outdir, "ssaha2Ref" )
        #system("mkdir -p %s" %refpath)
        #refpath = os.path.join( refpath, "ssaha2Ref")
        #system( "ssaha2Build -rtype 454 -save %s %s" %(refpath , options.ref) )
        #454
        #self.addChildTarget( Ls454Paired(options, refpath) )
       
        #Run bwa for abi_solid reads:
        #refpath = os.path.join(options.outdir, "bwaColorRef")
        #system("mkdir -p %s" %refpath)
        #refpath = os.path.join(refpath, "bwaColorRef")
        #system("bwa index -c -p %s %s" )

        #Run bowtie for abi_solid reads:
        if not self.options.ignoreSolid:
            refpath = os.path.join( options.outdir, "bowtieRef" )
            system("mkdir -p %s" %refpath)
            refpath = os.path.join( refpath, "bowtieRef")
            system( "bowtie-build -C %s %s" %(options.ref, refpath))
            #abi_solid
            self.addChildTarget( SolidSingle(options, refpath) )
            self.addChildTarget( SolidPaired(options, refpath) )

        #Combine results:
        if not self.options.nostats:
            self.setFollowOnTarget( CombineResults(options) )
        #self.setFollowOnTarget( AggregateResults(options) )

class CombineResults(Target):
    def __init__(self, options):
        Target.__init__(self)
        #Target.__init__(self, time=0.0025)
        self.options = options

    def run(self):
        outdir = self.options.outdir
        exps = ['illumina', 'ls454']
        if not self.options.ignoreSolid:
            exps.append('abi_solid')
        types = ['single', 'paired']
        outfile = os.path.join(outdir, "summaryStats.txt")
        header = ""
        outfh = open(outfile, "w")
        sumcounts = []
        for exp in exps:
            for type in types:
                expPath = os.path.join(outdir, exp, type)
                if os.path.exists(expPath):
                    f = open(os.path.join(expPath, "stats", "stats.txt"), "r")
                    i = 0
                    lastline = ""
                    for line in f.readlines():
                        if i == 0:
                            i = 1
                            if header == "":
                                header = line.strip()
                                outfh.write("%s\n" %header)
                                for c in range( len(header.split('\t')) -1 ):
                                    sumcounts.append(0)
                        elif re.match('Total', line):
                            lastline = line.strip()

                    if lastline != "":
                        items = lastline.split('\t')
                        if len(items) - 1 != len(sumcounts):
                        #if len(items) - 2 != len(sumcounts):
                            sys.stderr.write("file %s does not have enough stat fields. %d, %d\n" %( os.path.join(expPath, "stats.txt"), len(items) -1, len(sumcounts) ))
                        else:
                            outfh.write("%s-%s\t%s\n" %(exp, type, '\t'.join(items[1:])))
                            #for j in range(1, len(items)):
                            for j in range(1, len(sumcounts) + 1):
                                sumcounts[j -1] += int( items[j].split()[0] )
                    f.close()
        
        outfh.write("Total\t%d" %sumcounts[0])
        for c in range( 1, len(sumcounts) ):
            #outfh.write("\t%d" %(sumcounts[c]))
            outfh.write("\t%d (%.2f)" %(sumcounts[c], sumcounts[c]*100.0/sumcounts[0]))
        outfh.write("\n")

class SolidSingle(Target):
    def __init__(self, options, refpath):
        Target.__init__(self)
        #Target.__init__(self, time=0.0025)
        self.options = options
        self.refpath = refpath
    
    def run(self):
        solidSinglePath = os.path.join(self.options.readdir, "abi_solid", "single")
        if os.path.exists( solidSinglePath ):
            solidSingleOutdir = os.path.join(self.options.outdir, "abi_solid", "single")
            system("mkdir -p %s" %solidSingleOutdir)
            statsdir = os.path.join(solidSingleOutdir, "stats")
            system("mkdir -p %s" % statsdir)

            logger.info("Setting up solidSingle runs\n")
            for file in os.listdir( solidSinglePath ):
                if os.path.isdir(file):
                    continue
                self.addChildTarget( RunSolidSingle(self.options, solidSingleOutdir, self.refpath, solidSinglePath, file) )
            
            #Merge stat files:
            self.setFollowOnTarget( MergeResults(statsdir) )
            #Merge alignment files:
            #self.setFollowOnTarget( MergeBams(solidSingleOutdir) )

class SolidPaired(Target):
    def __init__(self, options, refpath):
        Target.__init__(self)
        #Target.__init__(self, time = 0.0025)
        self.options = options
        self.refpath = refpath
    
    def run(self):
        solidPairedPath = os.path.join(self.options.readdir, "abi_solid", "paired")
        if os.path.exists( solidPairedPath ):
            solidPairedOutdir = os.path.join(self.options.outdir, "abi_solid", "paired")
            system("mkdir -p %s" %solidPairedOutdir)
            statsdir = os.path.join(solidPairedOutdir, "stats")
            system("mkdir -p %s" %statsdir)

            logger.info("Setting up solidPaired runs\n")
            pairNames, nametail = getPairNames( os.listdir(solidPairedPath) )
            for pair in pairNames:
                if self.options.bamFormat:
                    self.addChildTarget( RunSolidPaired(self.options, solidPairedOutdir, self.refpath, solidPairedPath, "%s%s" %(pair, nametail), "%s%s" %(pair, nametail)) )
                else:
                    self.addChildTarget( RunSolidPaired(self.options, solidPairedOutdir, self.refpath, solidPairedPath, "%s_1%s" %(pair, nametail), "%s_2%s" %(pair, nametail)) )

            #Merge stat files:
            self.setFollowOnTarget( MergeResults(statsdir) )
            #Merge alignment files:
            #self.setFollowOnTarget( MergeBams(solidPairedOutdir) )

class Single(Target):
    def __init__(self, options, refpath, machine):
        Target.__init__(self)
        #Target.__init__(self, time = 0.0025)
        self.options = options
        self.refpath = refpath
        self.machine = machine
    
    def run(self):
        inPath = os.path.join(self.options.readdir, self.machine, "single")
        if os.path.exists( inPath ):
            outdir = os.path.join(self.options.outdir, self.machine, "single")
            system("mkdir -p %s" %outdir)
            statsdir = os.path.join(outdir, "stats")
            system("mkdir -p %s" %statsdir)

            logger.info("Setting up %sSingle runs\n" %self.machine)
            for file in os.listdir( inPath ):
                if os.path.isdir(file):
                    continue
                self.addChildTarget( RunBwaSingle(self.options, outdir, self.refpath, inPath, file) )
            
            #Merge stat files:
            self.setFollowOnTarget( MergeResults(statsdir) )
            #Merge alignment files:
            #self.setFollowOnTarget( MergeBams(outdir) )

class Paired(Target):
    def __init__(self, options, refpath, machine, maxInsert):
        Target.__init__(self)
        #Target.__init__(self, time=0.0025)
        self.options = options
        self.refpath = refpath
        self.machine = machine
        self.maxInsert = maxInsert
    
    def run(self):
        inPath = os.path.join(self.options.readdir, self.machine, "paired")
        if os.path.exists( inPath ):
            outdir = os.path.join(self.options.outdir, self.machine, "paired")
            system("mkdir -p %s" % outdir)
            statsdir = os.path.join(outdir, "stats")
            system("mkdir -p %s" %statsdir)

            logger.info("Setting up %sPaired runs\n" %self.machine)
            pairNames, nametail = getPairNames( os.listdir(inPath) )
            for pair in pairNames:
                if self.options.bamFormat:
                    self.addChildTarget( RunBwaPaired(self.options, outdir, self.refpath, inPath, "%s%s" %(pair, nametail), "%s%s" %(pair, nametail), self.maxInsert) )
                else:
                    self.addChildTarget( RunBwaPaired(self.options, outdir, self.refpath, inPath, "%s_1%s" %(pair, nametail), "%s_2%s" %(pair, nametail), self.maxInsert) )

            #Merge stat files:
            self.setFollowOnTarget( MergeResults(statsdir) )
            #Merge alignment files:
            #self.setFollowOnTarget( MergeBams(ls454PairedOutdir) )

class RunBwaPaired(Target):
    """
    """
    def __init__(self, options, outdir, refpath, indir, reads1, reads2, maxInsert):
        #Target.__init__(self, time=120)
        Target.__init__(self)
        self.outdir = outdir
        self.ref = refpath
        self.options = options
        self.indir = indir
        #self.name = reads1.split('_')[0]
        self.name = reads1.split('.')[0].split('_')[0]
        self.reads1 = reads1
        self.reads2 = reads2
        self.maxInsert = maxInsert

    def run(self):
        #HACK:
        if os.path.exists( os.path.join(self.outdir, "%s.bam" %self.name) ):
            return
        if os.path.exists( os.path.join(self.outdir, "inrangeOrUnmapped", "%s.bam" %self.name) ):
            return
        #END HACK


        localTempDir = self.getLocalTempDir()
        globalTempDir = self.getGlobalTempDir()
        #Copy ref to localTempDir:
        refpath = os.path.join(localTempDir, "ref")
        system("mkdir -p %s" %refpath)
        system("scp -C %s* %s" %(self.ref, refpath))
        #system("cp  %s* %s" %(self.ref, refpath))
        reffile = os.path.join(refpath, os.path.basename(self.ref) )

        #Copy input reads to localTempDir: 
        reads1 = os.path.join(localTempDir, self.reads1)
        system("scp -C %s %s" %( os.path.join(self.indir, self.reads1), reads1 ))
        #system("cp  %s %s" %( os.path.join(self.indir, self.reads1), reads1 ))
        reads2 = os.path.join(localTempDir, self.reads2)
        if self.reads1 != self.reads2:
            system("scp -C %s %s" %( os.path.join(self.indir, self.reads2), reads2 ))
            #system("cp %s %s" %( os.path.join(self.indir, self.reads2), reads2 ))

        sai1 = os.path.join(localTempDir, "1.sai")
        sai2 = os.path.join(localTempDir, "2.sai")
        
        if self.options.bamFormat:
            system("bwa aln -t %d %s -b1 %s > %s" %(self.options.numBwaThreads, reffile, reads1, sai1) )
            system("bwa aln -t %d %s -b2 %s > %s" %(self.options.numBwaThreads, reffile, reads2, sai2) )
        else:
            system("bwa aln -t %d %s %s > %s" %(self.options.numBwaThreads, reffile, reads1, sai1) )
            system("bwa aln -t %d %s %s > %s" %(self.options.numBwaThreads, reffile, reads2, sai2) )
        
        bamfile = os.path.join(localTempDir, "%s-Mapped.bam" %self.name)
        #system("bwa sampe -n 100 -a %d %s %s %s %s %s > %s" % (self.maxInsert, reffile, sai1, sai2, reads1, reads2, samfile))
        system("bwa sampe -n 100 -a %d %s %s %s %s %s | samtools view -S -b - > %s" % (self.maxInsert, reffile, sai1, sai2, reads1, reads2, bamfile))
            
        globalBamfile = os.path.join(globalTempDir, "%s.bam" %self.name)
        system("mv %s %s" %(bamfile, globalBamfile))

        self.addChildTarget( RunStats(self.options, self.name, "paired", self.outdir, globalTempDir) )
        
        #Cleanup:
        #system("rm -R %s" %localTempDir)
        #system("rm -R %s" %globalTempDir)

class RunBwaSingle(Target):
    """
    """
    def __init__(self, options, outdir, refpath, indir, reads):
        #Target.__init__(self, time=120)
        Target.__init__(self)
        self.outdir = outdir
        self.ref = refpath
        self.options = options
        self.indir = indir
        self.name = reads.split('.')[0]
        self.reads = reads

    def run(self):
        #HACK:
        if os.path.exists( os.path.join(self.outdir, "%s.bam" %self.name) ):
            return
        if os.path.exists( os.path.join(self.outdir, "inrangeOrUnmapped", "%s.bam" %self.name) ):
            return
        #END HACK
        
        localTempDir = self.getLocalTempDir()
        globalTempDir = self.getGlobalTempDir()
        #Copy ref to localTempDir:
        refpath = os.path.join(localTempDir, "ref")
        system("mkdir -p %s" %refpath)
        system("mkdir -p %s" %refpath)
        system("cp %s* %s" %(self.ref, refpath))
        reffile = os.path.join(refpath, os.path.basename(self.ref) )

        #Copy input reads to localTempDir: 
        reads = os.path.join(localTempDir, self.reads)
        system("scp -C %s %s" %( os.path.join(self.indir, self.reads), reads ))
        #system("cp %s %s" %( os.path.join(self.indir, self.reads), reads ))

        sai = os.path.join(localTempDir, "%s.sai" %self.name)
        if self.options.bamFormat:
            system("bwa aln -t %d %s -b %s > %s" %(self.options.numBwaThreads, reffile, reads, sai) )
        else:
            system("bwa aln -t %d %s %s > %s" %(self.options.numBwaThreads, reffile, reads, sai) )

        bamfile = os.path.join(localTempDir, "%s-Mapped.bam" %self.name)
        system("bwa samse -n 100 %s %s %s | samtools view -S -b - > %s " % ( reffile, sai, reads, bamfile ))
        globalBamfile = os.path.join(globalTempDir, "%s.bam" %self.name)
        system("mv %s %s" %(bamfile, globalBamfile))

        self.addChildTarget( RunStats(self.options, self.name, "single", self.outdir, globalTempDir) )
        
        #Cleanup:
        #system("rm -R %s" %localTempDir)
        #system("rm -R %s" %globalTempDir)

class RunSolidSingle(Target):
    """
    """
    def __init__(self, options, outdir, refpath, indir, reads):
        #Target.__init__(self, time=120)
        Target.__init__(self)
        self.outdir = outdir
        self.ref = refpath
        self.options = options
        self.indir = indir
        self.name = reads.split('.')[0]
        self.reads = reads

    def run(self):
        #HACK:
        if os.path.exists( os.path.join(self.outdir, "%s.bam" %self.name) ):
            return
        if os.path.exists( os.path.join(self.outdir, "inrangeOrUnmapped", "%s.bam" %self.name) ):
            return
        #END HACK
        
        localTempDir = self.getLocalTempDir()
        globalTempDir = self.getGlobalTempDir()
        #Copy ref to localTempDir:
        refpath = os.path.join(localTempDir, "ref")
        system("mkdir -p %s" %refpath)
        system("cp %s* %s" %(self.ref, refpath))
        reffile = os.path.join(refpath, os.path.basename(self.ref) )

        #Copy input reads to localTempDir: 
        reads = os.path.join(localTempDir, self.reads)
        system("scp -C %s %s" %( os.path.join(self.indir, self.reads), reads ))
        #system("cp %s %s" %( os.path.join(self.indir, self.reads), reads ))

        bamfile = os.path.join(localTempDir, "%s-Mapped.bam" %self.name)
        system("bowtie -p 3 -C -S -t %s %s | samtools view -S -b - > %s" %(reffile, reads, bamfile)) 
        globalBamfile = os.path.join(globalTempDir, "%s.bam" %self.name)
        system("mv %s %s" %(bamfile, globalBamfile))
        
        self.addChildTarget( RunStats(self.options, self.name, "single", self.outdir, globalTempDir) )
        
        #Cleanup:
        #system("rm -R %s" %localTempDir)
        #system("rm -R %s" %globalTempDir)

class RunSolidPaired(Target):
    """
    """
    def __init__(self, options, outdir, refpath, indir, reads1, reads2):
        #Target.__init__(self, time=120)
        Target.__init__(self)
        self.outdir = outdir
        self.ref = refpath
        self.options = options
        self.indir = indir
        self.name = reads1.split('_')[0]
        self.reads1 = reads1
        self.reads2 = reads2

    def run(self):
        #HACK:
        if os.path.exists( os.path.join(self.outdir, "%s.bam" %self.name) ):
            return
        if os.path.exists( os.path.join(self.outdir, "inrangeOrUnmapped", "%s.bam" %self.name) ):
            return
        #END HACK
        
        localTempDir = self.getLocalTempDir()
        globalTempDir = self.getGlobalTempDir()
        #Copy ref to localTempDir:
        refpath = os.path.join(localTempDir, "ref")
        system("mkdir -p %s" %refpath)
        system("cp %s* %s" %(self.ref, refpath))
        reffile = os.path.join(refpath, os.path.basename(self.ref))

        #Copy input reads to localTempDir: 
        reads1 = os.path.join(localTempDir, self.reads1)
        #system("cp %s %s" %( os.path.join(self.indir, self.reads1), reads1 ))
        system("scp -C %s %s" %( os.path.join(self.indir, self.reads1), reads1 ))
        reads2 = os.path.join(localTempDir, self.reads2)
        #system("cp %s %s" %( os.path.join(self.indir, self.reads2), reads2 ))
        system("scp -C %s %s" %( os.path.join(self.indir, self.reads2), reads2 ))

        bamfile = os.path.join(localTempDir, "%s-Mapped.bam" %self.name)
        system("bowtie -p 3 -C -S -t %s -1 %s -2 %s |samtools view -S -b - > %s" %(reffile, reads1, reads2, bamfile)) 
        globalBamfile = os.path.join(globalTempDir, "%s.bam" %self.name)
        system("mv %s %s" %(bamfile, globalBamfile))

        self.addChildTarget( RunStats(self.options, self.name, "paired", self.outdir, globalTempDir) )
        
        #Cleanup:
        #system("rm -R %s" %localTempDir)
        #system("rm -R %s" %globalTempDir)

class RunStats(Target):
    def __init__( self, options, name, libLayout, outdir, globalTempDir ):
        #Target.__init__(self, time=120)
        Target.__init__(self)
        self.name = name
        self.options = options
        self.libLayout = libLayout
        self.outdir = outdir
        self.indir = globalTempDir

    def run(self):
        globalTempDir = self.indir
        globalBamfile = os.path.join(globalTempDir, "%s.bam" %self.name)
        localTempDir = self.getLocalTempDir()
        bamfile = os.path.join(localTempDir, "%s.bam" %self.name)
        system("cp %s %s" %(globalBamfile, bamfile))
        
        #getStats:
        sortedbamFile = os.path.join(localTempDir, "%s-sorted.bam" %self.name)
        if self.options.nostats:
            sortPrefix = os.path.join(localTempDir, "%s-sorted" % self.name)
            system("samtools sort -n %s %s" %(bamfile, sortPrefix))
        else: 
            sortedbamFile, statsFile = stats(localTempDir, self.name, bamfile)
            statsdir = os.path.join(self.outdir, "stats")
            system("cp %s %s" %(statsFile, statsdir))

        if self.options.getMappedReads2:
            inrangeDir = os.path.join(self.outdir, "inrangeOrUnmapped")
            system("mkdir -p %s" %(inrangeDir))
            outfilePrefix = os.path.join(localTempDir, "%s" %self.name)
            filterReads2(localTempDir, sortedbamFile, outfilePrefix, self.options.readsRange)
            #system("cp %s.bam %s" %(outfilePrefix, inrangeDir))
            system("scp -C %s.bam %s" %(outfilePrefix, inrangeDir))
        
        if self.options.getMappedReads:
            #Filter out readPairs that have both reads unmapped:
            outfilePrefix = os.path.join(localTempDir, "%s" %self.name)
            if self.libLayout == "single":
                filterReadsSingle(sortedbamFile, outfilePrefix)
            else:
                filterReads(localTempDir, sortedbamFile, outfilePrefix)
            #Move filtered alignment file back to output directory
            #system("cp %s.bam %s" %(outfilePrefix,self.outdir))
            system("scp -C %s.bam %s" %(outfilePrefix,self.outdir))

class MergeResults(Target):
    def __init__(self, indir):
        #Target.__init__(self, time=0.25)
        Target.__init__(self)
        self.indir = indir

    def run(self):
        pattern = ".+-stats.txt"
        files = getfiles(pattern, self.indir)
        if len(files) == 0:
            return
        outfile = os.path.join(self.indir, "stats.txt")
        outfh = open(outfile, "w")
        outfh.write("Name\tTotal\tFailure\tDuplicates\tMapped\tPairedInSequencing\tRead1\tRead2\tProperlyPaired\tWithItselfAndMateMapped\tSingletons\tMateMappedToDiffChr\tMateMappedToDiffChr(mapQ>=5)\t\
                     UniquelyMapped\tWithItselfAndMateUniquelyMapped\tWithItselftAndMateUniquelyMappedAndProperlyPaired\tUniquelyMappedMateNotUniquelyMapped\tWithItselfAndMateNotUniquelyMapped\n")
        totalcounts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for file in files:
            counts  = []
            f = open(file, 'r')
            for line in f.readlines():
                items = line.strip().split()
                if len(items) < 2:
                    continue
                counts.append( int(items[0]) )
            f.close()
            total = counts[0]
            name = os.path.basename(file).rstrip('-stats.txt')
            outfh.write("%s\t%d" %(name, total))
            totalcounts[0] += total
            for i in range(1,len(counts)):
                totalcounts[i] += counts[i]
                outfh.write( "\t%d (%.2f%%)" %(counts[i], counts[i]*100.0/total) )
            outfh.write("\n")

        outfh.write("Total\t%d" %totalcounts[0])
        for i in range(1,len(totalcounts)):
            outfh.write( "\t%d (%.2f%%)" %(totalcounts[i], totalcounts[i]*100.0/totalcounts[0]) )
        outfh.write("\n")
        outfh.close()

class MergeBams(Target):
    def __init__(self, indir):
        Target.__init__(self)
        self.indir = indir

    def run(self):
        files = os.listdir(self.indir)
        if len(files) == 0:
            return
        filesStr = " ".join([ os.path.join(self.indir, file) for file in files ])
	logger.info("Merging %s\n" %filesStr)
        mergeFile = os.path.join(self.indir, "merge.bam")
        mergeCmd = "samtools merge %s %s" %(mergeFile, filesStr)
        system( mergeCmd )
        sortPrefix = os.path.join(self.indir, "sortedMerge")
	logger.info("Sorted %s with prefix %s\n" %(mergeFile, sortPrefix))
        sortCmd = "samtools sort -n %s %s" %( mergeFile, sortPrefix)
        system(sortCmd)
        #Convert to sam:
        system("samtools view %s.bam -o %s.sam" %(sortPrefix, sortPrefix))

def getfiles(pattern, indir):
    files = []
    #print "GETFILES %s\n" %indir
    if not os.path.isdir(indir):
        print "%s is NOT DIR\n" % indir
    
    filelist = os.listdir( indir )
    filelist.sort()
    for file in filelist:
        if re.search( pattern, file ):
            files.append( os.path.join( indir, file ) )
    return files

def filterReads(indir, file, outfilePrefix):
    #reffile = os.path.join(indir, "ref.fa")
    #system("cp %s %s" %(ref, reffile))

    mapMap = os.path.join(indir, "mm.bam")
    system("samtools view -F12 -b -o %s %s" %(mapMap, file) )
    
    mapUnmap = os.path.join(indir, "mu.bam")
    system("samtools view -f8 -F4 -b -o %s %s" %(mapUnmap, file) )

    unmapMap = os.path.join(indir, "um.bam")
    system("samtools view -f4 -F8 -b -o %s %s" %(unmapMap, file) )

    #Merge to outfile:
    mergeFile = os.path.join(indir, "merge.bam")
    system("samtools merge %s %s %s %s" %(mergeFile, mapMap, mapUnmap, unmapMap) )
    system( "rm %s %s %s" %(mapMap, mapUnmap, unmapMap) )
    #Sort by name:
    system("samtools sort -n %s %s" %(mergeFile, outfilePrefix))
    return

def filterReadsSingle(file, outfilePrefix):
    system("samtools view -F4 -b %s | samtools sort -n - %s" %(file, outfilePrefix))
    return

#====== Filter to only get reads that mapped to the input range and unmapped reads ======
def filterReads2(indir, file, outfilePrefix, posrange):
    #sort by position:
    prefix = os.path.join(indir, "sorted")
    system("samtools sort %s %s" %(file, prefix))
    system("samtools index %s.bam" % prefix)

    inrange = os.path.join(indir, "inrange.bam")
    system("samtools view -F4 -b %s.bam %s > %s" %(prefix, posrange, inrange))
    inrangeSort = os.path.join(indir, "inrangeSorted")
    system("samtools sort -n %s %s" %(inrange, inrangeSort))

    unmapped = os.path.join(indir, "unmapped.bam")
    system("samtools view -f4 -b %s > %s" %(file, unmapped))

    #Merge to outfile:
    mergeFile = os.path.join(indir, "merge.bam")
    system("samtools merge %s %s.bam %s" %(mergeFile, inrangeSort, unmapped))
    system("rm %s %s %s.bam" %(inrange, unmapped, inrangeSort))

    #Sort by name:
    system("samtools sort -n %s %s" %(mergeFile, outfilePrefix))
    return

def getPairNames( files ):
    names = []
    tail = ""
    for file in files:
        if os.path.isdir(file):
            continue
        items = file.split("_")
        if items[0] not in names:
            names.append( items[0] )
            extra = '_'.join(items[1:])
            extra = extra.lstrip('12')
            if tail == "":
                tail = extra
            elif tail != extra:
                sys.stderr.write('Input files must have the same name format! Found different ones: *%s *%s\n' %(tail, extra))
                sys.exit(1)
    return names, tail

#def decodeFlag(flag):
#    flagDict = {"isPaired": 0, "properlyPaired":1, "unmapped":2, \
#                "unmappedMate":3, "strand":4, "mateStrand":5, "firstRead":6,\
#                "secondRead":7, "notPrimaryAlign":8, "failedQual":9,\
#                "duplicate":10}
#
#    for key in flagDict:
#        if 2**flagDict[key] & flag == 0:
#            flagDict[key] = False
#        else:
#            flagDict[key] = True
#    return flagDict

#def getStats(file, outfile):
#    f = open(file, 'r')
#    pair = 0
#    prevReadUniquelyMapped = False
#    prevFlagDict = {}
#    counts = {"uniquely_mapped":0,\
#              "with_itself_and_mate_uniquely_mapped":0,\
#              "with_itself_and_mate_uniquely_mapped_and_properly_paired":0,\
#              "uniquely_mapped_mate_not_uniquely_mapped":0,\
#              "with_itself_and_mate_not_uniquely_mapped":0}
#
#    #for line in f.readlines():
#    #while (line = f.readline()) != '':
#    for line in f:
#        pair += 1
#        line = line.strip()
#        if len(line) == 0 or line[0] == '@':
#            continue
#        items = line.split('\t')
#        flagDict = decodeFlag( int(items[1]) )
#        uniquelyMapped = True
#        if len(items) >= 12 and re.search( 'XA:', '\t'.join([i for i in items[11:]]) ):
#            uniquelyMapped = False
#        
#        if uniquelyMapped:
#            if not flagDict['unmapped']:#Current read is mapped uniquely
#                counts["uniquely_mapped"] += 1
#                
#                if flagDict['isPaired'] and pair %2 == 0:
#                    if prevReadUniquelyMapped:
#                        counts["with_itself_and_mate_uniquely_mapped"] += 2
#                        if flagDict['properlyPaired']:
#                            counts["with_itself_and_mate_uniquely_mapped_and_properly_paired"] += 2
#                    elif not prevFlagDict['unmapped']:#prev read mapped but not uniquely
#                        counts["uniquely_mapped_mate_not_uniquely_mapped"] += 1
#        else:
#            if not flagDict["unmapped"]:#current read is mapped but Not uniquely
#                if flagDict["isPaired"] and pair %2 == 0:
#                    if prevReadUniquelyMapped:
#                        counts["uniquely_mapped_mate_not_uniquely_mapped"] += 1
#                    elif not prevFlagDict['unmapped']:#prev read is mapped but not uniquely:
#                        counts["with_itself_and_mate_not_uniquely_mapped"] +=2
#
#        prevReadUniquelyMapped = uniquelyMapped and (not flagDict["unmapped"])
#        prevFlagDict = flagDict
#    f.close()
#
#    #Print out stats:
#    outfh = open(outfile, 'w')
#    for k in sorted( counts.keys() ):
#        outfh.write("%d\t%s\n" %(counts[k], k))
#    outfh.close()

def stats(indir, name, file):
    #Convert sam to bam for flagstat
    #bamfile = os.path.join(indir, "%s.bam" %name)
    #system("samtools view -b -o %s -S %s" %(bamfile, file))
    
    #Sort bamfile by name
    sortPrefix = os.path.join(indir, "%s-sorted" %name)
    system("samtools sort -n %s %s" %(file, sortPrefix))
    
    flagstatFile = os.path.join(indir, "%s-flagstat.txt" %name)
    system("samtools flagstat %s.bam > %s" %(sortPrefix, flagstatFile))

    #Convert sorted bam back to sam (sorted sam) to parse and find uniquely-mapped statistics
    #system("samtools view %s.bam -o %s.sam" %(sortPrefix, sortPrefix))
    uniqueStatFile = os.path.join(indir, "%s-ustat.txt" %name)
    #getStats("%s.sam" %sortPrefix, uniqueStatFile)
    system("parseBam.py < %s.bam > %s" %(sortPrefix, uniqueStatFile))

    statsFile = os.path.join(indir, "%s-stats.txt" %name)
    #MergeStats:
    system("cat %s %s > %s" %(flagstatFile, uniqueStatFile, statsFile))
    #Return the sorted samfile, and statsfile
    return "%s.bam" % sortPrefix, statsFile

def initOptions( parser ):
    #parser.add_option('-bwape', dest='bwape', help='Comma separated list of input bam/fastq read-files: reads1.bam,reads2.bam,otherReads1.bam,otherReads2.bam. Reads2 must follow reads1 file immediately in the list.')
    #parser.add_option('-bwa', dest='bwa', help='Comma separated list of input bam/fastq files.')
    #parser.add_option('-bwasw' dest='bwasw', help='Comma separated list of input bam/fastq files.')
    #parser.add_option('-ssaha2', dest='ssaha2', help='Comma separated list of input fastq files.')
    #parser.add_option('-ssaha2pe', dest='ssaha2pe', help='Comma separated list of input fastq files.')
    parser.add_option('-o', '--outdir', dest='outdir', default = ".", help='Output directory. Default is current directory')
    parser.add_option('-i', '--readdir', dest='readdir', help='Required argument, directory where input fastq files are. The directory should have this hierachy: readdir/sampleName/sequencingPlatform/pairedOrSingle/file*.fastq')
    #parser.add_option('--illuminaInsertSize', dest='iInsertSize', default = 200, type ="int", help="Insert size of illumina paired reads. Default = 200")
    parser.add_option('--illuminaMaxInsertSize', dest='iMaxInsertSize', default = 1000, type = "int", help="Maximum allowed insert size of illumina paired reads. Default = 1000")
    parser.add_option('--454MaxInsertSize', dest='ls454MaxInsertSize', default = 5000, type = "int", help="Maximum allowed insert size of illumina paired reads. Default = 5000")
    #parser.add_option('--454InsertSizeMin', dest='ls454InsertSizeMin', default=1750, type="int", help="Insert size of 454 paired reads. Default=1750")
    #parser.add_option('--454InsertSizeMax', dest='ls454InsertSizeMax', default=3750, type="int", help="Insert size of 454 paired reads. Default=3750")
    parser.add_option('--numBwaThreads', dest='numBwaThreads', type = "int", default = 1, help = "Default =1")
    parser.add_option('--getMappedReads', dest='getMappedReads', action='store_true', default=False, help="If specified, get reads that mapped or have mapped mates")
    parser.add_option('--getMappedReads2', dest='getMappedReads2', action='store_true', default=False, help="Ignore this")
    parser.add_option('--nostats', dest='nostats', action='store_true', default=False, help='If specified, produce no mapping statistics')
    parser.add_option('--readsRange', dest='readsRange', default='chr6:28,477,797-33,428,773', help="Ignore this")
    parser.add_option('--ignoreSolid', dest='ignoreSolid', action='store_true', default=False)
    parser.add_option('--bam', dest='bamFormat', action='store_true', default=False, help="Must specify if input files are in bam format. Otherwise it is assumed that input files are in fastq format")
    parser.add_option('--bwaRefIndex', dest='bwaRefIndex', help="Specify bweRef prefix if have already done the indexing for the reference")

def checkOptions( options, args, parser ):
    if len(args) == 0 and options.bwaRefIndex == None:
        parser.error("No reference sequence or reference index was given.\n")
    elif len(args) > 0 and not os.path.exists(args[0]):
        parser.error("Reference file %s does not exist.\n" %args[0])
    elif len(args) > 0:
        options.ref = args[0]
    
    if not os.path.exists( options.outdir ):
        system("mkdir -p %s" %options.outdir)

    if options.readdir == None:
        parser.error("Please specify readdir\n")

    #if options.bwape:
    #    options.bwape = options.bwape.split(',')
    #if options.bwa:
    #    options.bwa = options.bwa.split(',')
    #if options.bwasw:
    #    options.bwasw = options.bwasw.split(',')
    #if options.ssaha2:
    #    options.ssaha2 = options.ssaha2.split(',')
    #if options.ssaha2pe:
    #    options.ssaha2pe = options.ssaha2pe.split(',')

def main():
    usage = "Usage: %prog [options] ref.fa"
    parser = OptionParser( usage = usage )
    Stack.addJobTreeOptions(parser)

    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( options, args, parser )

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

    
if __name__ == "__main__":
    from referenceViz.src.mapReadsToRef import *
    main()
