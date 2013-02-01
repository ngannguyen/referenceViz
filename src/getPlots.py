#!/usr/bin/env python

"""
nknguyen at soe dot ucsc dot edu
May 31 2011
Script to generate different analyses plots for different experiments (different cactus parameters) 
of the reference project.
Analyses included are: contiguity,coverage,n50,snp,indeldist,indeltab,cnv,snpcheck,indelcheck,snpcheckplot,indelcheckplot,indel,
all
"""

import os, sys, re, time
import xml.etree.ElementTree as ET
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import setLogLevel
from sonLib.bioio import getTempDirectory

class Setup(Target):
    """Set up analysis runs for all experiments
    """
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.indir = options.indir
        self.outdir = options.outdir
        self.analyses = options.analyses
        self.pdflatex = options.pdflatex
        self.dbsnp = options.dbsnp
        self.pgsnp = options.pgsnp
        self.dbindel = options.dbindel
        self.pgindel = options.pgindel
        self.refstart = options.refstart
        self.refend = options.refend
        self.indelMaxSize = options.indelMaxSize
        self.filter= options.filter
        self.pileupSnp = options.pileupSnp

    def run(self):
         experiments = os.listdir(self.indir)
         system("mkdir -p %s" %self.outdir)
         for exp in experiments:
             filteredSamples = getFilteredSamples(exp)

             indir = os.path.join(self.indir, exp)
             if not os.path.isdir(indir):
                 continue

             outdir = os.path.join(self.outdir, exp)
             system("mkdir -p %s" %outdir)
             
             if 'contiguity' in self.analyses:
                 pattern1 = "contiguityStats_.+\.xml"
                 self.addChildTarget( Contiguity(indir, outdir, pattern1, filteredSamples) )
                 self.addChildTarget( Contiguity(indir, outdir, pattern1, filteredSamples, includeCoverage=True) )
                 pattern2 = "contiguityStatsNoDuplication_.+\.xml"
                 self.addChildTarget( Contiguity(indir, outdir, pattern2, filteredSamples, includeCoverage=True) )
             
             if 'coverage' in self.analyses:
                 self.addChildTarget( Coverage(indir, outdir, filteredSamples) )
             
             if 'n50' in self.analyses:
                 pattern = "pathStats_[^_]+\.xml"
                 self.addChildTarget( N50(indir, outdir, filteredSamples, pattern) )
                 pattern = "pathStats_ignoreAdjacencies_[^_]+\.xml"
                 self.addChildTarget( N50(indir, outdir, filteredSamples, pattern) )
             
             if 'snp' in self.analyses:
                 #Snp plots
                 pattern = "snpStats_[^_]+\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern, filteredSamples) )
                 pattern = "snpStats_filtered_[^_]+\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern, filteredSamples) )
                 pattern = "snpStats_[^_]+_recurrent\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern, filteredSamples) )
                 pattern = "snpStats_filtered_[^_]+_recurrent\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern, filteredSamples) )
             
             if 'snpcheck' in self.analyses:
                 pattern = "snpStats_hg19\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_hg19_withAggregates\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_hg19_recurrent\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_hg19_recurrent_withAggregates\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_filtered_hg19\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_filtered_hg19_withAggregates\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_filtered_hg19_recurrent\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
                 pattern = "snpStats_filtered_hg19_recurrent_withAggregates\.xml"
                 self.addChildTarget( SnpCheck(indir, outdir, pattern, filteredSamples, self.dbsnp, self.pgsnp, self.refstart, self.refend, self.filter, self.pileupSnp) )
             
             if 'indelcheck' in self.analyses:
                 pattern = "pathStats_hg19.xml"
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 0, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 5, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 pattern = "pathStats_hg19_withAggregates.xml"
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 0, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 5, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 pattern = "pathStats_ignoreAdjacencies_hg19.xml"
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 0, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 5, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 pattern = "pathStats_ignoreAdjacencies_hg19_withAggregates.xml"
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 0, self.indelMaxSize, self.refstart, self.refend, self.filter) )
                 self.addChildTarget( IndelCheck(indir, outdir, pattern, filteredSamples, self.dbindel, self.pgindel, 5, self.indelMaxSize, self.refstart, self.refend, self.filter) )
             
             if 'indeldist' in self.analyses:
                 self.addChildTarget( IndelDist(indir, outdir, filteredSamples) )

             if 'indel' in self.analyses:
                 self.addChildTarget( Indel(indir, outdir, filteredSamples) )

             if 'indeltab' in self.analyses:
                 pathPatterns = [ "pathStats_[^_]+\.xml", "pathStats_ignoreAdjacencies_[^_]+\.xml" ]
                 snpPatterns = ["snpStats_[^_]+\.xml", "snpStats_[^_]+_recurrent\.xml", "snpStats_filtered_[^_]+\.xml", "snpStats_filtered_[^_]+_recurrent\.xml"]
                 for p1 in pathPatterns:
                     for p2 in snpPatterns:
                         self.addChildTarget( IndelTab(indir, outdir, self.pdflatex, filteredSamples, p1, p2) )
             if 'cnv' in self.analyses:
                 self.addChildTarget( Cnv(indir, outdir, filteredSamples) )
             
         self.setFollowOnTarget( CheckPlot(self.analyses, self.outdir) )

class CheckPlot(Target):
    def __init__(self, analyses, outdir):
        Target.__init__(self, time=0.25)
        self.outdir = outdir
        self.analyses = analyses
    
    def run(self):
         experiments = os.listdir(self.outdir)
         for exp in experiments:
             outdir = os.path.join(self.outdir, exp)
             if not os.path.isdir(outdir):
                 continue
             
             if 'indelcheckplot' in self.analyses:
                 self.addChildTarget( IndelCheckPlot(outdir, outdir) )
             
             if 'snpcheckplot' in self.analyses:
                 self.addChildTarget( SnpCheckPlot(outdir, outdir) )

class Contiguity(Target):
    """
    """
    def __init__(self, indir, outdir, pattern, filteredSamples, includeCoverage=False):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern
        self.includeCoverage = includeCoverage
        self.filteredSamples = filteredSamples

    def run(self):
        #pattern = "contiguityStats_.+\.xml"
        files = getfiles(self.pattern, self.indir)
        filesStr = " ".join(files)
       
        if len(files) >= 1:
            cmd = "contiguityPlot.py %s --outdir %s " %(filesStr, self.outdir) 
            if self.includeCoverage:
                cmd += " --includeCoverage "
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class Coverage(Target):
    """
    """
    def __init__(self, indir, outdir, filteredSamples):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.filteredSamples = filteredSamples

    def run(self):
        infile = os.path.join(self.indir, "coverageStats.xml")
        if os.path.exists( infile ):
            cmd = "coveragePlot.py %s --outdir %s " %(infile, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class N50(Target):
    def __init__(self, indir, outdir, filteredSamples, pattern):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.filteredSamples = filteredSamples
        self.pattern = pattern

    def run(self):
        #pattern = "pathStats_[^_]+\.xml"
        files = getfiles(self.pattern, self.indir)
        prefix = 'n50'
        if len(files) > 0:
            prefix += os.path.basename( files[0] ).lstrip('pathStats_').rstrip('.xml').rstrip('\.hg19')

        filesStr = " ".join(files)
        sortkey = "scaffoldPathN50"
        keys = "sequenceN50,blockN50,contigPathN50,scaffoldPathN50"
        #keys = "blockN50,contigPathN50,scaffoldPathN50"
        if len(files) >=1:
            cmd = "n50Plot.py %s --sortkey %s --keys %s --outdir %s -p %s" %(filesStr, sortkey, keys, self.outdir, prefix)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class Snp(Target):
    def __init__(self, indir, outdir, pattern, filteredSamples):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern
        self.filteredSamples = filteredSamples

    def run(self):
        #pattern = "snpStats_.+\.xml"
        #pattern = "snpStatsIntersection_.+\.xml"
        files = getfiles(self.pattern, self.indir)
        filesStr = " ".join(files)
        if len(files) >=1:
            cmd = "snpPlot.py %s --outdir %s" %(filesStr, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class SnpCheckPlot(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        system("snpCheckPlot.py --indir %s --outdir %s" %(self.indir, self.outdir))

class IndelCheckPlot(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        system("indelCheckPlot.py --indir %s --outdir %s" %(self.indir, self.outdir))

class SnpCheck(Target):
    def __init__(self, indir, outdir, pattern, filteredSamples, dbsnp, pgsnp, refstart, refend, filter, pileupSnp):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern
        self.filteredSamples = filteredSamples
        self.dbsnp = dbsnp
        self.pgsnp = pgsnp
        self.refstart = refstart
        self.refend = refend
        self.filter = filter
        self.pileupSnp = pileupSnp

    def run(self):
        files = getfiles(self.pattern, self.indir)
        filters = {None:''}
        if self.filter != None:
            filters[self.filter] = '_noRepeats'
        
        for file in files:
            name = os.path.basename(file).lstrip('snpStats_').rstrip('\.xml')
            for filter in filters:
                outfile = os.path.join(self.outdir, "dbsnpCheck%s_%s.txt" %(filters[filter], name))
                if self.pgsnp == None:
                    if self.refstart != None and self.refend != None:
                        if filter == None:
                            system("snpStats.py -s %d -e %d %s %s > %s" %(self.refstart, self.refend, file, self.dbsnp, outfile))
                        else:
                            system("snpStats.py -f %s -s %d -e %d %s %s > %s" %(self.filter, self.refstart, self.refend, file, self.dbsnp, outfile))
                    else:
                        if filter == None:
                            system("snpStats.py %s %s > %s" %(file, self.dbsnp, outfile))
                        else:
                            system("snpStats.py -f %s %s %s > %s" %(self.filter, file, self.dbsnp, outfile))
                else:
                    if self.refstart != None and self.refend != None:
                        if filter == None:
                            system("snpStats.py --pileupSnp %s -s %d -e %d --pgSnp %s %s %s > %s" %(self.pileupSnp, self.refstart, self.refend, self.pgsnp, file, self.dbsnp, outfile))
                        else:
                            system("snpStats.py --pileupSnp %s -f %s -s %d -e %d --pgSnp %s %s %s > %s" %(self.pileupSnp, self.filter, self.refstart, self.refend, self.pgsnp, file, self.dbsnp, outfile))
                    else:
                        if filter == None:
                            system("snpStats.py --pileupSnp %s --pgSnp %s %s %s > %s" %(self.pileupSnp, self.pgsnp, file, self.dbsnp, outfile))
                        else:
                            system("snpStats.py --pileupSnp %s -f %s --pgSnp %s %s %s > %s" %(self.pileupSnp, self.filter, self.pgsnp, file, self.dbsnp, outfile))
                
class IndelCheck(Target):
    def __init__(self, indir, outdir, pattern, filteredSamples, dbsnp, pgsnp, wobble, cutoff, refstart, refend, filter):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern
        self.filteredSamples = filteredSamples
        self.dbsnp = dbsnp
        self.pgsnp = pgsnp
        self.wobble = wobble
        self.cutoff = cutoff
        self.refstart = refstart
        self.refend = refend
        self.filter = filter

    def run(self):
        files = getfiles(self.pattern, self.indir)
        filters = {None:''}
        if self.filter != None:
            filters[self.filter] = '_noRepeats'

        for file in files:
            name = os.path.basename(file).lstrip('pathStats_').rstrip('\.xml')
            for filter in filters.keys():
                outfile = os.path.join( self.outdir, "indelCheck%s_%d_%s.txt" %(filters[filter], self.wobble, name) )
                if self.pgsnp == None:
                    if self.refstart != None and self.refend != None:
                        if filter == None:
                            system("indelStats.py -s %d -e %d -c %d -w %d %s %s > %s" %(self.refstart, self.refend, self.cutoff, self.wobble, file, self.dbsnp, outfile))
                        else:
                            system("indelStats.py -f %s -s %d -e %d -c %d -w %d %s %s > %s" %(self.filter, self.refstart, self.refend, self.cutoff, self.wobble, file, self.dbsnp, outfile))
                    else:
                        if filter == None:
                            system("indelStats.py -c %d -w %d %s %s > %s" %(self.cutoff, self.wobble, file, self.dbsnp, outfile))
                        else:
                            system("indelStats.py -f %s -c %d -w %d %s %s > %s" %(self.filter, self.cutoff, self.wobble, file, self.dbsnp, outfile))
                else:
                    if self.refstart != None and self.refend != None:
                        if filter == None:
                            system("indelStats.py -s %d -e %d -c %d -w %d --pgSnp %s %s %s > %s" %(self.refstart, self.refend, self.cutoff, self.wobble, self.pgsnp, file, self.dbsnp, outfile))
                        else:
                            system("indelStats.py -f %s -s %d -e %d -c %d -w %d --pgSnp %s %s %s > %s" %(self.filter, self.refstart, self.refend, self.cutoff, self.wobble, self.pgsnp, file, self.dbsnp, outfile))
                    else:
                        if filter == None:
                            system("indelStats.py -c %d -w %d --pgSnp %s %s %s > %s" %(self.cutoff, self.wobble, self.pgsnp, file, self.dbsnp, outfile))
                        else:
                            system("indelStats.py -f %s -c %d -w %d --pgSnp %s %s %s > %s" %(self.filter, self.cutoff, self.wobble, self.pgsnp, file, self.dbsnp, outfile))

class Indel(Target):
    def __init__(self, indir, outdir, filteredSamples):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.filteredSamples = filteredSamples

    def run(self):
        pattern = "pathStats_[^_]+_withAggregates.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)

        if len(files) >=1:
            cmd = "indelPlot.py %s --outdir %s " %(filesStr, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class IndelDist(Target):
    def __init__(self, indir, outdir, filteredSamples):
        Target.__init__(self, time = 0.25)
        self.indir = indir
        self.outdir = outdir
        self.filteredSamples = filteredSamples
    
    def run(self):
        pattern = "pathStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)

        if len(files) >=1:
            cmd = "indelDistPlot.py %s --outdir %s " %(filesStr, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

class IndelTab(Target):
    def __init__(self, indir, outdir, pdflatex, filteredSamples, pathstatsPattern, snpstatsPattern):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.pdflatex = pdflatex
        self.filteredSamples = filteredSamples
        self.pathPattern = pathstatsPattern
        self.snpPattern = snpstatsPattern

    def run(self):
        #pattern = "pathStats_[^_]+\.xml"
        files1 = getfiles(self.pathPattern, self.indir)
        filesStr = " ".join(files1)
        #pattern = "snpStats_[^_]+\.xml"
        files2 = getfiles(self.snpPattern, self.indir)
        filesStr += " " + " ".join(files2)

        outPrefix = "indelStats"
        if len(files1) > 0 and len(files2) > 0:
            basename1 = os.path.basename( files1[0] ).lstrip('pathStats_').rstrip('\.xml').rstrip('hg19').strip('reference')
            basename2 = os.path.basename( files2[0] ).lstrip('snpStats_').rstrip('\.xml').strip('hg19').strip('reference')
            outPrefix = "indelStats_%s_%s" %( basename1, basename2 )

        if len(files1) >=1 and len(files2) >=1:
            cmd = "indelTable.py %s --outPrefix %s --outdir %s " %(filesStr, outPrefix, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)
            #Get make pdf for tex files:
            #if self.pdflatex:
            #    prefix = "indelStats_"
            #    outfile = os.path.join( self.outdir, prefix + "hg19-reference" + ".tex" )
            #    system( "pdflatex %s" %(outfile) )
            #    system( "mv %s %s" %("indelStats_*.pdf", self.outdir) )
            #    system( "rm -f *.aux *.log")
                #for file in files:
                    #m = re.match( "pathStats_(.+)\.xml", os.path.basename(file) )
                    #outfile = os.path.join( self.outdir, prefix + m.group(1) + ".tex" )

class Cnv(Target):
    def __init__(self, indir, outdir, filteredSamples):
        Target.__init__(self, time=0.25)
        self.indir = indir
        self.outdir = outdir
        self.filteredSamples = filteredSamples

    def run(self):
        infile = os.path.join(self.indir, "copyNumberStats.xml")
        if os.path.exists( infile ):
            cmd = "cnvPlot.py %s --outdir %s " %(infile, self.outdir)
            if self.filteredSamples != "":
                cmd += " --filteredSamples %s" %(self.filteredSamples)
            system(cmd)

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
             
def getFilteredSamples(exp):
    #filtered sample:
    filteredSamples = ''
    items = exp.split('_')
    if len(items) > 2:
        filteredSamples = items[ len(items) -1 ]
    return filteredSamples

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Required. Location of all the experiments') 
    parser.add_option('-o', '--outdir', dest='outdir', help='Required. Output directory')
    parser.add_option('-l', '--pdflatex', dest='pdflatex', action="store_true", default=False, help='If specified, convert tex files to pdf using "pdflatex" ("pdflatex" must be installed)')
    parser.add_option('--dbsnp', dest='dbsnp', help='dbSnps file') 
    parser.add_option('--pgsnp', dest='pgsnp', help='pgSnps file') 
    parser.add_option('--dbindel', dest='dbindel', help='dbSnps indels file') 
    parser.add_option('--pgindel', dest='pgindel', help='pgSnps indels file') 
    parser.add_option('-a', '--analyses', dest='analyses', default='all', 
                      help='Comma separated string of different analyses to perform.'
                      'Analyses are within the list:[contiguity,coverage,n50,snp,indeldist,indeltab,cnv,snpcheck,indelcheck,snpcheckplot,indelcheckplot,indel,all].'
                      'The default string is "all", which means all the analyses included.'
                      'contiguity = Contiguity; '
                      'coverage: Coverage; '
                      'n50: N50; '
                      'snp: SNP rate; '
                      'indel, indeltab: Indel rate; '
                      'indeldist: Indel length distribution; '
                      'snpcheck, snpcheckplot, indelcheck, indelcheckplot: dbSNP and 1000 Genomes project SNPs and short indels (<= 10bp) validation; '
                      'cnv: copy number variation')
    parser.add_option('--indelMaxSize', dest='indelMaxSize', type='int', default=10, help="Only indels with size <= than this cutoff are included in the comparisons with dbSnps indels. Default = 10")
    parser.add_option('-r', '--ref', dest='ref', help='current reference (GRCh37/hg19) sequence')
    parser.add_option('-f', '--filter', dest='filter', help='File contain regions to ignore in the stats (format:chr\\tchromStart\\tchromEnd). Will ignore all snps lie within this region. Default=no filtering')
    parser.add_option('--pileupSnp', dest='pileupSnp')

def checkOptions( args, options, parser ):
    if not options.indir:
        parser.error('Location of input experiments is required and not given.\n')
    if not options.outdir:
        parser.error('Output direcotry is required but not given.\n')
    
    if re.search('all', options.analyses):
        options.analyses = 'contiguity,coverage,n50,snp,indeldist,indeltab,cnv,snpcheck,indelcheck,indelcheckplot,snpcheckplot,indel'
    options.analyses = (options.analyses).split(',')
    alist = ['contiguity', 'coverage', 'n50', 'snp', 'indeldist', 'indeltab', 'cnv', 'snpcheck', 'indelcheck', 'indelcheckplot', 'snpcheckplot', 'indel']
    for a in options.analyses:
        if a not in alist:
            parser.error('Analysis %s is unknown\n' % a)
    
    if 'snpcheck' in options.analyses and (not options.dbsnp or not os.path.exists(options.dbsnp)):
        parser.error('Dbsnp file does not exist or was not provided\n')
    
    if options.pgsnp and not os.path.exists(options.pgsnp):
        parser.error('pgsnp file does not exist\n')
    
    if 'indelcheck' in options.analyses and (not options.dbindel or not os.path.exists(options.dbindel)):
        parser.error('Dbindel file does not exist or was not provided\n')
    if options.pgindel and not os.path.exists(options.pgindel):
        parser.error('pgsnp indel file does not exist\n')

    options.refstart = None
    options.refend = None
    if options.ref != None:
        f = open(options.ref, 'r')
        for l in f:
            if len(l) >0 and l[0] == '>':
                items = l.split('.')
                if len(items) >=6 and items[5] == '1':
                    s = int(items[3])
                    e = s + int(items[4])
                    if options.refstart == None or s < options.refstart:
                        options.refstart = s
                    if options.refend == None or options.refend < e:
                        options.refend = e
        f.close()
    

def main():
    usage = ( 'usage: %prog [options]\n'
              '%prog is a pipeline to generate various analyses plots for different experiments.\n'
              'Analyses included are: contiguity, coverage, n50, SNPs, indels, indel distribution, cnv, dbSNP validation (SNP and short indels)') 

    parser = OptionParser( usage = usage )
    Stack.addJobTreeOptions(parser)

    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    
    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)


if __name__ == "__main__":
    from referenceViz.src.getPlots import *
    main()

