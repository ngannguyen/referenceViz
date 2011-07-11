#!/usr/bin/env python

"""
Create Snp plots
nknguyen at soe dot ucsc dot edu
Jun 15 2011
"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

#from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.font_manager import FontProperties


class Sample():
    def __init__( self, xmlnode ):
        self.name = xmlnode.attrib[ 'sampleName' ]
        self.refname = xmlnode.attrib[ 'referenceName' ]
        self.totalErrors = int( xmlnode.attrib[ 'totalErrors' ] )
        self.totalCalls = int( xmlnode.attrib[ 'totalCalls' ] )
        self.errPerSite = float(self.totalErrors)/self.totalCalls
        #self.totalSites = int( xmlnode.attrib[ 'totalSites' ] )
        #self.totalCorrect = int( xmlnode.attrib[ 'totalCorrect' ] )

def readfiles( options ):
    statsList = []
    for f in options.files:
        samples = []
        xmltree = ET.parse( f )
        root = xmltree.getroot()
        for s in root.findall( 'statsForSample' ):
            name = s.attrib[ 'sampleName' ]
            if name != 'ROOT' and name != '':
                samples.append( Sample( s ) )
        statsList.append( samples )
    return statsList

def initOptions( parser ):
    parser.add_option( '--outdir', dest='outdir', default='.', help='Output directory' ) 

def checkOptions( args, options, parser ):
    if len( args ) < 1:
        parser.error( 'Please provide at least one snpStats xml file.\n' )
    options.files = []
    for f in args:
        if not os.path.exists(f):
            parser.error( 'File %s does not exist.\n'  % f )
        options.files.append( f )

def drawSnpData( axes, samples ):
    ydata = []
    xlabels = []
    for s in samples:
        ydata.append( s.errPerSite )
        xlabels.append( s.name )
    axes.plot( ydata, marker='.', markersize=10.0, linestyle='none' )
    libplot.editSpine( axes )
    pyplot.xlabel( 'Samples' )
    pyplot.ylabel( 'Snps per site' )

    return xlabels

def drawSnpPlot( options, samples ):
    samples = sorted( samples, key=lambda s:s.errPerSite, reverse=True )
    if len( samples ) < 1:
        return
    refname = samples[0].refname
    options.out = os.path.join( options.outdir, 'snp_' + refname )
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    xlabels = drawSnpData( axes, samples )
    title = 'SNPs'
    axes.set_title( title )

    #set ticks
    axes.set_xticks( range( 0, len(xlabels) ) )
    axes.set_xticklabels( xlabels )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation( 90 )

    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    axes.set_xlim( -0.5, len(samples) -0.5 )
    
    libplot.writeImage( fig, pdf, options )
   
def getSample( samples, name ):
    for s in samples:
        if s.name == name:
            return s
    return None

def drawCompareSnpData( axes, xsamples, ysamples ):
    xdata = []
    ydata = []
    for x in xsamples:
        y = getSample( ysamples, x.name )
        if y != None:
            xdata.append( x.errPerSite )
            ydata.append( y.errPerSite )
    
    axes.plot( xdata, ydata, marker='.', markersize=10.0, linestyle='none' )
    
    xmax = max( xdata )
    ymax = max( ydata )
    maxval = max( [xmax, ymax] ) 

    #Draw y=x line:
    xl = [-maxval*0.1, maxval*1.1]
    yl = [-maxval*0.1, maxval*1.1]
    axes.plot(xl, yl, color="0.9")

    libplot.editSpine( axes )
    pyplot.xlabel( xsamples[0].refname )
    pyplot.ylabel( ysamples[0].refname )

    return maxval

def drawCompareSnpPlot( options, xsamples,  ysamples ):
    if len(xsamples) < 1 or len(ysamples) < 1:
        return
    xrefname = xsamples[0].refname
    yrefname = ysamples[0].refname
    options.out = os.path.join( options.outdir, 'snp_' + xrefname + '_' + yrefname )
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    maxval = drawCompareSnpData( axes, xsamples, ysamples )
    title = 'Snps %s vs %s' % (xrefname, yrefname)
    axes.set_title( title )

    libplot.setTicks( axes )
    axes.set_xlim( -maxval*0.1, maxval*1.1 )
    axes.set_ylim( -maxval*0.1, maxval*1.1 )
    libplot.writeImage( fig, pdf, options )

def main():
    usage = ('Usage: %prog [options] file1.xml file2.xml\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    statsList = readfiles( options )

    for samples in statsList:
        drawSnpPlot( options, samples )

    if len(statsList) >= 2:
        for i in range( len(statsList) -1 ):
            for j in range( i+1, len(statsList) ):
                drawCompareSnpPlot( options, statsList[i], statsList[j] )

    

if __name__ == "__main__":
    main()
