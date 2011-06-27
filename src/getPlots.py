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

    def run(self):
         experiments = os.listdir(self.indir)
         system("mkdir -p %s" %self.outdir)
         for exp in experiments:
             indir = os.path.join(self.indir, exp)
             if not os.path.isdir(indir):
                 continue

             outdir = os.path.join(self.outdir, exp)
             system("mkdir -p %s" %outdir)
             self.addChildTarget( Contiguity(indir, outdir) )
             self.addChildTarget( Coverage(indir, outdir) )
             self.addChildTarget( N50(indir, outdir) )
             self.addChildTarget( Snp(indir, outdir) )
             self.addChildTarget( IndelDist(indir, outdir) )
             self.addChildTarget( IndelTab(indir, outdir) )
        #Cleanup:
        #self.setFollowOnTarget( Cleanup(self.output) )

class Contiguity(Target):
    """
    """
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        #if self.indir == 'mhcHumanVariantsNsRemovedV0/.gitignore':
        #    return
        pattern = "contiguityStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)
       
        if len(files) >= 1:
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
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        pattern = "snpStats_.+\.xml"
        files = getfiles(pattern, self.indir)
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
    def __init__(self, indir, outdir):
        Target.__init__(self, time=0.00025)
        self.indir = indir
        self.outdir = outdir

    def run(self):
        pattern = "pathStats_.+\.xml"
        files = getfiles(pattern, self.indir)
        filesStr = " ".join(files)
        if len(files) >=1:
            system("indelTable.py %s --outdir %s" %(filesStr, self.outdir))

def getfiles(pattern, indir):
    files = []
    print "GETFILES %s\n" %indir
    if not os.path.isdir(indir):
        print "%s is NOT DIR\n" % indir
        
    for file in os.listdir( indir ):
        if re.search( pattern, file ):
            files.append( os.path.join( indir, file ) )
    return files

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Required. Location of all the experiments') 
    parser.add_option('-o', '--outdir', dest='outdir', help='Required. Output directory')

def checkOptions( args, options, parser ):
    if not options.indir:
        parser.error('Location of input experiments is required and not given.\n')
    if not options.outdir:
        parser.error('Output direcotry is required but not given.\n')

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

