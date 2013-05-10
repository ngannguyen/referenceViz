#!/usr/bin/env python

"""
Tue Apr 23 15:49:19 PDT 2013
Compute SNPs between each pair of samples for all samples.
   and SNPs bewtween each sample to the Reference (e.g C.Ref)
Create Snp plots

Input: sampleList, referenceName, expPath, outputdir

nknguyen at soe dot ucsc dot edu
"""
import os, sys, re
from optparse import OptionParser

import cPickle as pickle
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system
from sonLib.bioio import logger

from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.ktserverLauncher import KtserverLauncher
from cactus.progressive.multiCactusProject import MultiCactusProject

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('PDF')
import libPlotting as libplot
import matplotlib.pyplot as pyplot
import matplotlib.pylab as pylab
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties


#######################################################
#=================== OBJECTS =========================#
#######################################################
class Sample():
    def __init__( self, xmlnode ):
        self.name = xmlnode.attrib[ 'sampleName' ]
        self.refname = xmlnode.attrib[ 'referenceName' ]
        self.totalErrors = int( xmlnode.attrib[ 'totalErrors' ] )
        self.totalCalls = int( xmlnode.attrib[ 'totalCalls' ] )
        self.errPerSite = 0
        if self.totalCalls != 0:
            self.errPerSite = float(self.totalErrors)/self.totalCalls
        #self.totalSites = int( xmlnode.attrib[ 'totalSites' ] )
        #self.totalCorrect = int( xmlnode.attrib[ 'totalCorrect' ] )

#######################################################
#=================== PIPELINE ========================#
#######################################################
class Setup(Target):
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        if self.options.snpdir:
            self.setFollowOnTarget( ProcessSnpStats(self.options, self.options.snpdir) )
        else:
            #launch ktserver
            if self.options.launchKtserver:
                exp = getExperiment(self.options.expPath)
                self.options.ktserver = launchDbServer(exp)

            #run snpStats with each sample as reference for all samples
            snpdir = os.path.join(self.options.outdir, "snpStats")
            system("mkdir -p %s" %snpdir)
            samples = self.options.samples + [self.options.ref]
            for i in xrange(0, len(samples)):
                refsam = samples[i]
                snpfile = os.path.join(snpdir, "snpStats_%s.xml" %(refsam))
                self.addChildTarget( RunSnpStats(self.options.expPath, refsam, snpfile) )
            
            #after got all the snpStats necessary, process the data and make plots
            self.setFollowOnTarget( ProcessSnpStats(self.options, snpdir) )

class RunSnpStats(Target):
    '''Run the snpStats script for the input ref
    '''
    def __init__(self, expPath, ref, outfile):
        Target.__init__(self)
        self.expPath = os.path.abspath(expPath)
        self.ref = ref
        self.outfile = outfile

    def run(self):
        exp = getExperiment(self.expPath)
        dbstr = exp.getDiskDatabaseString()
        #ktserver = launchDbServer(exp)
        cmd = "snpStats --cactusDisk '%s' --outputFile %s --referenceEventString %s" %(dbstr, self.outfile, self.ref)
        system(cmd)
        #killDbServer(ktserver, exp)

class ProcessSnpStats(Target):
    '''read all snpStats of all the references
    '''
    def __init__(self, options, snpdir):
        Target.__init__(self)
        self.options = options
        self.snpdir = snpdir

    def run(self):
        if not self.options.snpdir and self.options.launchKtserver: #kill ktserver:
            exp = getExperiment(self.options.expPath)
            killDbServer(self.options.ktserver, exp)
        ref2samples = readfiles( self.options.samples, self.snpdir )
        drawSnpPlot(self.options, ref2samples)
        
#######################################################
#============ UTILITIES FUNCTIONS ====================#
#######################################################
def drawSnpPlot(options, ref2samples):
    #All the samples sorted indecreasing order of SNP rate (relative to the reference)
    refname = options.ref
    samples_ref = sorted( ref2samples[refname].values(), key=lambda s:s.errPerSite, reverse=True )
    sampleNames = [s.name for s in samples_ref]
    if len(sampleNames) < 1:
        return
    
    #Get SNP data w.r.t the reference
    refYdata = [ s.errPerSite for s in samples_ref ]
    refXdata = range(1, len(refYdata) + 1)
    logger.info( "minRef %d\n" % min(refYdata) )
    logger.info( "maxRef %d\n" % max(refYdata) )

    #Get average of SNP data of each sample w.r.t each of the other samples
    ydataList = []
    for name in sampleNames:
        data = []
        for ref in sampleNames:
            if ref == name:
                continue
            if ref not in ref2samples or name not in ref2samples[ref]:
                sys.stderr.write("Ref %s not in ref2samples or name %s not in ref2samles[ref]\n" %(ref, name))
                sys.exit(1)
            sample = ref2samples[ref][name]
            data.append( sample.errPerSite )
        ydataList.append(data)

    #Calculate the average and std:
    meanYdata = [np.mean(data) for data in ydataList]
    stdYdata = [np.std(data) for data in ydataList]
    offset = 0.25
    xdata = [x + offset for x in refXdata]

    logger.info( "minMean %d\n" % min(meanYdata) )
    logger.info( "maxMean %d\n" % max(meanYdata) )

    #Set up
    options.out = os.path.join(options.outdir, "snpPlot")
    fig, pdf = libplot.initImage( 8.0, 12.0, options )
    axes = fig.add_axes([0.15, 0.08, 0.78, 0.85])

    #Plot: #switch x and y axes because there are so many samples
    linenames = [refname, "other"]
    colors =["#1F78B4", "#E31A1C", "#4DAF4A"] 
    #SNPs w.r.t the reference
    l1 = axes.plot( refYdata, refXdata, color=colors[0], marker='.', markersize=12.0, linestyle='none')
    #SNPs w.r.t other samples:
    l2 = axes.errorbar( meanYdata, xdata, xerr=stdYdata, color=colors[1], markeredgecolor=colors[1], markersize=12.0, fmt='.')

    #Draw horizontal lines for each samples
    upperYdata = [y + stdYdata[i] for i, y in enumerate(meanYdata)]
    lowerYdata = [y - stdYdata[i] for i, y in enumerate(meanYdata)]
    miny = min( [min(refYdata), min(lowerYdata)] )
    maxy = max( [max(refYdata), max(upperYdata)] )
    rangey = maxy - miny
    maxy += rangey*0.02
    miny -= rangey*0.02
    for x in xdata:
        x -= offset/2
        axes.plot([miny, maxy], [x, x], color="#515151", linestyle='-', linewidth=0.25 )
    axes.set_xlim( miny, maxy )
    axes.set_ylim( 0.5, len(sampleNames) + 0.5 )

    #Legend
    libplot.editSpine(axes)
    fontP = FontProperties()
    fontP.set_size("x-small")
    legend = axes.legend([l1, l2], linenames, 'best', numpoints=1, prop=fontP)
    legend._drawFrame = False

    sampleLabels = [libplot.properName(l) for l in sampleNames]
    pyplot.yticks([x - offset/2 for x in xdata], sampleLabels, fontproperties=fontP)
    axes.set_ylabel( 'Samples' )
    axes.set_xlabel( 'SNPs per site' )
    title = 'SNPs'
    axes.set_title( title )
    
    libplot.writeImage( fig, pdf, options )
   
def readfiles( samples, indir ):
    ref2samples = {}
    for f in os.listdir( indir ):
        ref = f.split('.')[0].split('_')[-1]
        ref2samples[ref] = {}#key = sampleName, val = Sample

        file = os.path.join(indir, f)
        xmltree = ET.parse( file )
        root = xmltree.getroot()
        for s in root.findall( 'statsForSample' ):
            name = s.attrib[ 'sampleName' ]
            #if name != 'aggregate' and name != 'ROOT' and name != '' and 
            if name in samples:
                ref2samples[ref][name] =  Sample(s)
    return ref2samples

def readSampleList(file):
    samples = []
    f = open(file, 'r')
    for line in f:
        samples.append(line.strip())
    f.close()
    return samples

def getExperiment(expPath):
    expXml = ET.parse(expPath).getroot()
    return ExperimentWrapper(expXml)

def launchDbServer(experiment):
    ktserver = None
    if experiment.getDbType() == "kyoto_tycoon":
        ktserver = KtserverLauncher()
        ktserver.spawnServer(experiment, readOnly=True)
    if ktserver:
        logger.info("Got ktserver, not NULL")
    else:
        logger.info("Did not launch ktserver")
    return ktserver

def killDbServer(ktserver, experiment):
    if experiment.getDbType() == "kyoto_tycoon":
        ktserver.killServer(experiment)

def initOptions( parser ):
    #Input: sampleList, referenceName, projectPath, outputdir
    parser.add_option( '-e', '--experimentPath', dest='expPath', help='Require argument. Experiment path. Default=%default' )
    parser.add_option( '-s', '--samples', dest='samples', help='Required argument. File containing list of samples. Default=%default' )
    parser.add_option( '-r', '--reference', dest='ref', help='Required argument. Reference sequence Name. Default=%default')
    parser.add_option( '-o', '--outdir', dest='outdir', default='.', help='Output directory' ) 
    parser.add_option( '--snpStatsDir', dest='snpdir', help='If specified, assumed that all the snpStats have already been computed and are located in this directory')
    parser.add_option("--ktserver", dest="launchKtserver", action='store_true', default=False, help='If ktserver has not been launched, this option must be specified')
    #parser.add_option( '--numOutliners', dest='numOutliners', default=1, help='Number of outliners' ) 
    #parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    #parser.add_option('--outgroup', dest='outgroup', help='Name of outgroup sample')

def checkOptions( args, options, parser ):
    if not options.expPath or not options.samples or not options.ref:
        parser.error("Please specify all required arguments: expPath, samples and ref\n")
    if not os.path.exists( options.expPath ):
        parser.error("ProjectPath %s does not exist\n" %options.expPath)
    options.samples = readSampleList( options.samples )

def main():
    usage = ('Usage: %prog [options]\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from referenceViz.src.snpPlot_all2all import *
    main()
