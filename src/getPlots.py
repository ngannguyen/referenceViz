#!/usr/bin/env python

"""
nknguyen at soe dot ucsc dot edu
May 31 2011
Script to generate different analysis plots for different experiments (different cactus parameters) 
of the reference project.
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

    def run(self):
         experiments = os.listdir(self.indir)
         system("mkdir -p %s" %self.outdir)
         for exp in experiments:
             indir = os.path.join(self.indir, exp)
             if not os.path.isdir(indir):
                 continue

             outdir = os.path.join(self.outdir, exp)
             system("mkdir -p %s" %outdir)
             if 'contiguity' in self.analyses:
                 pattern1 = "contiguityStats_.+\.xml"
                 self.addChildTarget( Contiguity(indir, outdir, pattern1) )
                 self.addChildTarget( Contiguity(indir, outdir, pattern1, includeCoverage=True) )
                 pattern2 = "contiguityStatsNoDuplication_.+\.xml"
                 self.addChildTarget( Contiguity(indir, outdir, pattern2, includeCoverage=True) )
             if 'coverage' in self.analyses:
                 self.addChildTarget( Coverage(indir, outdir) )
             if 'n50' in self.analyses:
                 self.addChildTarget( N50(indir, outdir) )
             if 'snp' in self.analyses:
                 pattern = "snpStats_.+\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern) )
                 pattern = "snpStatsIntersection_.+\.xml"
                 self.addChildTarget( Snp(indir, outdir, pattern) )
             if 'indeldist' in self.analyses:
                 self.addChildTarget( IndelDist(indir, outdir) )
             if 'indeltab' in self.analyses:
                 self.addChildTarget( IndelTab(indir, outdir, self.pdflatex) )
             if 'cnv' in self.analyses:
                 self.addChildTarget( Cnv(indir, outdir) )
        #Cleanup:
        #self.setFollowOnTarget( Cleanup(self.output) )

class Contiguity(Target):
    """
    """
    def __init__(self, indir, outdir, pattern, includeCoverage=False):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern
        self.includeCoverage = includeCoverage

    def run(self):
        #pattern = "contiguityStats_.+\.xml"
        files = getfiles(self.pattern, self.indir)
        filesStr = " ".join(files)
       
        if len(files) >= 1:
            if self.includeCoverage:
                system("contiguityPlot.py %s --outdir %s --includeCoverage" %(filesStr, self.outdir) )
            else:
                system("contiguityPlot.py %s --outdir %s" %(filesStr, self.outdir) )


class Coverage(Target):
    """
    """
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        infile = os.path.join(self.indir, "coverageStats.xml")
        if os.path.exists( infile ):
            system("coveragePlot.py %s --outdir %s" %(infile, self.outdir))

class N50(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        pattern = "pathStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)
        sortkey = "scaffoldPathN50"
        keys = "sequenceN50,blockN50,contigPathN50,scaffoldPathN50"
        #keys = "blockN50,contigPathN50,scaffoldPathN50"
        if len(files) >=1:
            system("n50Plot.py %s --sortkey %s --keys %s --outdir %s" %(filesStr, sortkey, keys, self.outdir))

class Snp(Target):
    def __init__(self, indir, outdir, pattern):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.pattern = pattern

    def run(self):
        #pattern = "snpStats_.+\.xml"
        #pattern = "snpStatsIntersection_.+\.xml"
        files = getfiles(self.pattern, self.indir)
        filesStr = " ".join(files)
        if len(files) >=1:
            system("snpPlot.py %s --outdir %s" %(filesStr, self.outdir))

class IndelDist(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self, time = 0.00025)
        self.indir = indir
        self.outdir = outdir
    
    def run(self):
        pattern = "pathStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)

        if len(files) >=1:
            system("indelDistPlot.py %s --outdir %s" %(filesStr, self.outdir))

class IndelTab(Target):
    def __init__(self, indir, outdir, pdflatex):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir
        self.pdflatex = pdflatex

    def run(self):
        pattern = "pathStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)
        pattern = "snpStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr += " " + " ".join(files)

        if len(files) >=1:
            system("indelTable.py %s --outdir %s" %(filesStr, self.outdir))
            #Get make pdf for tex files:
            if self.pdflatex:
                prefix = "indelStats_"
                outfile = os.path.join( self.outdir, prefix + "hg19-reference" + ".tex" )
                system( "pdflatex %s" %(outfile) )
                system( "mv %s %s" %("indelStats_*.pdf", self.outdir) )
                system( "rm -f *.aux *.log")
                #for file in files:
                    #m = re.match( "pathStats_(.+)\.xml", os.path.basename(file) )
                    #outfile = os.path.join( self.outdir, prefix + m.group(1) + ".tex" )

class Cnv(Target):
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        infile = os.path.join(self.indir, "copyNumberStats.xml")
        if os.path.exists( infile ):
            system("cnvPlot.py %s --outdir %s" %(infile, self.outdir))

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

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Required. Location of all the experiments') 
    parser.add_option('-o', '--outdir', dest='outdir', help='Required. Output directory')
    parser.add_option('-l', '--pdflatex', dest='pdflatex', action="store_true", default=False, help='If specified, convert tex files to pdf using "pdflatex" (This must be installed)')
    parser.add_option('-a', '--analyses', dest='analyses', default='all', \
                      help='Comma separated string of different analyses to perform.\n\
                      Analyses are within the list:[contiguity,coverage,n50,snp,indeldist,indeltab,cnv,all].\n\
                      The default string is "all", which means all the analyses included.')

def checkOptions( args, options, parser ):
    if not options.indir:
        parser.error('Location of input experiments is required and not given.\n')
    if not options.outdir:
        parser.error('Output direcotry is required but not given.\n')
    if re.search('all', options.analyses):
        options.analyses = 'contiguity,coverage,n50,snp,indeldist,indeltab,cnv'
    options.analyses = (options.analyses).split(',')
    alist = ['contiguity', 'coverage', 'n50', 'snp', 'indeldist', 'indeltab', 'cnv']
    
    for a in options.analyses:
        if a not in alist:
            parser.error('Analysis %s is unknown\n' % a)


def main():
    usage = ( 'usage: %prog [options]'
              '%prog is a pipeline to generate various analysis plots for different (inputed) experiments') 

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
    #_test()
    main()

