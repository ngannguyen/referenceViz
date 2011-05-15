#!/usr/bin/env python

"""
Create contiguity plots
nknguyen at soe dot ucsc dot edu
May 11 2011
Input: two contiguityStats.xml files to be compared and contrast
(e.g: contiguityStats_reference.xml and contiguityStats_hg19.xml)
"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator

class Bucket:
    def __init__( self, bucketElement ):
        self.start = int( bucketElement.attrib[ 'from' ] )
        self.end = int( bucketElement.attrib[ 'to' ] )
        self.mid = ( self.start + self.end )/2
        self.correct = int( bucketElement.attrib[ 'correct' ] )
        self.samples = int( bucketElement.attrib[ 'samples' ] )
        self.correctPerSample = float( bucketElement.attrib[ 'correctPerSample' ] )
        self.cumulativeCorrect = int( bucketElement.attrib[ 'cumulativeCorrect' ] )
        self.cumulativeSamples = int( bucketElement.attrib[ 'cumulativeSamples' ] )
        self.cumulativeCorrectPerSample = float( bucketElement.attrib[ 'cumulativeCorrectPerSample' ] )
        
class Sample( list ):
    def __init__( self, sampleElement ):
        self.name = sampleElement.attrib[ 'sampleName' ]
        self.reference = sampleElement.attrib[ 'referenceName' ]
        for bucket in sampleElement.findall( 'bucket' ):
            self.append( Bucket( bucket ) )

class Stats( list ): #each Stats represents one input XML file
    def __init__( self, name ):
        self.name = name
    def setRefName( self, refname ):
        self.refname = refname

def getSample( stats, name ):
    for sample in stats:
        if sample.name == name:
            return sample
    return None

def getColors( num ):
    dimsize = num ** (1/3.0) 
    if int(dimsize) < dimsize:
        dimsize = int(dimsize) + 1
    else:
        dimsize = int(dimsize)
    dimsize += 2 #don't use black (0,0,0) and white (1,1,1)
    #Return a list with 'num' number of colors
    colors = []  #each element is a tuple of 3 floats representing R, G, B values, Eg: (0.1, 0.5, 0.8)
    
    vals = linspace( 0, 1, dimsize )
    for r in vals:
        for g in vals:
            for b in vals:
                colors.append( (r, g, b) )

    return colors

def setTicks( axes ):
    axes.xaxis.set_ticks_position('bottom')
    axes.yaxis.set_ticks_position('left')
    minorLocator = LogLocator( base=10, subs = range(1,10) )
    axes.xaxis.set_minor_locator( minorLocator )

def setAxisLimits( axes, ycutoff ):
    axes.set_xscale('log')
    axes.set_ylim( ycutoff, 1.001 )

def drawLegend( lines, sampleNames, options ):
    if not options.legendElements:
        legend = pyplot.legend( lines, sampleNames)
        legend._drawFrame=False
    elif len(lines) == len(options.legendElements):
        legend = pyplot.legend( lines, options.legendElements )
        legend._drawFrame=False
    else:
        sys.stderr.write('Number of items in --legendElements is different '
                         'from the number of lines plotted\n' )
    
def editSpine( axes ):    
    for loc, spine in axes.spines.iteritems():
        if loc in ['left', 'bottom']:
            spine.set_position( ('outward', 10) )
        elif loc in ['right', 'top']:
            spine.set_color('none')
        else:
            raise ValueError( 'Unknown spine location %s\n' % loc )


def drawData( axes, stats ):
    halfsize = len(stats)/2 + len(stats)%2
    colors = getColors( halfsize )
    styles = { 0:'-', 1:'--' }

    dash = 0
    colorindex = 0
    lines = []
    sampleNames = []

    for sample in stats:
        sampleNames.append(sample.name)
        xdata = []
        ydata = []
        for bucket in sample:
            xdata.append( bucket.mid )
            ydata.append( bucket.correctPerSample )
        
        if not dash:
            colorindex += 1
        
        l = axes.plot( xdata, ydata, color=colors[colorindex], linestyle=styles[dash], linewidth=2.0 )
        lines.append(l)
        
        dash = not dash
    
    editSpine( axes )
    return lines, sampleNames

def drawCompareData( axes, xstats, ystats ):
    #Only draw the overlapped samples
    colors = getColors( len(xstats) )
    colorindex = 0
    lines = []
    sampleNames = []
    for xsample in xstats:
        if xsample == ystats.refname:
            continue
        ysample = getSample( ystats, xsample.name )
        if ysample is None:
            continue
        if len(xsample) != len(ysample):
            sys.stderr.write( "Error: Two xml files do not have the same number of buckets for sample %s\n" % xsample.name )
            sys.exit( 1 )
        xdata = []
        ydata = []
        colorindex += 1
        for i in range( len( xsample ) ): #each bucket
            if xsample[i].mid != ysample[i].mid:
                sys.stderr.write( "Two xml files have different buckets\n " )
                sys.exit( 1 )
            xdata.append( xsample[i].correctPerSample )
            ydata.append( ysample[i].correctPerSample )
        
        l = axes.plot( xdata, ydata, color=colors[colorindex], marker='.', markersize=9.0, linestyle='none' )
        lines.append( l )
        sampleNames.append( xsample.name )
        #pyplot.line.set_label( xsample.name )
    editSpine( axes )
    return lines, sampleNames


def setAxes( fig ):
    return fig.add_axes( [0.12, 0.1, 0.83, 0.85] )

def drawContiguityPlot( options, stats ):
    options.out = stats.refname #name of output file
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = setAxes( fig )
    
    lines, sampleNames = drawData( axes, stats )
    drawLegend( lines, sampleNames, options )
    setAxisLimits( axes, options.ycutoff )
    setTicks( axes )

    libplot.writeImage( fig, pdf, options )

def drawCompareContiguityPlot( options, xstats, ystats ):
    options.out = xstats.refname + "_" + ystats.refname
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = setAxes( fig )

    lines, sampleNames = drawCompareData( axes, xstats, ystats )
    legend = pyplot.legend( lines, sampleNames, numpoints = 1)
    legend._drawFrame=False
    #axes.set_xlim( 0, 1.001 )
    #axes.set_ylim( 0, 1.001 )
    #setAxisLimits( axes, options.ycutoff )
    setTicks( axes )
    
    libplot.writeImage( fig, pdf, options )

def readfiles( options ):
    statsList = [] #each element represents one input XML file (one contiguity plot)
    for f in options.files:
        name = os.path.basename( f ).split( '.' )[0]
        stats = Stats( name )

        xmltree = ET.parse( f )
        root = xmltree.getroot()
        for sample in root.findall( 'statsForSample' ):
            name = sample.attrib[ 'sampleName' ]
            if name != '' and name != 'ROOT':
                stats.append( Sample( sample ) )
        if len(stats) > 0:
            stats.setRefName( stats[0].reference )

        statsList.append( stats )

    return statsList


def initOptions( parser ):
    parser.add_option('--title', dest='title', default='Contiguous Statistics',
                       help='Based title of the plots, default=%default')
    parser.add_option('--legendElements', dest='legendElements',
                       help='Specify the legend text - comma separated list' )
    parser.add_option('--ycutoff', dest='ycutoff', default=0.5, type='float',
                       help='Only points with y-value from ycutoff to 1 are displayed')

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error('Please specify at least one contiguityStats file.\n')
    options.files = []
    for f in args:
        if not os.path.exists( f ):
            parser.error('%s does not exist\n' %f)
        options.files.append( os.path.abspath( f ) )

    if options.legendElements:
        options.legendElements = options.legendElements.split(',')

def main():
    usage = ( 'usage: %prog [options] file1.xml file2.xml\n\n'
              '%prog takes in contiguityStats.xml files and create an image file' )

    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()

    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    statsList = readfiles( options )

    for stats in statsList:
        drawContiguityPlot( options, stats )
        
    if len(statsList) >= 2:
        for i in range( len(statsList) -1 ):
            for j in range( i + 1, len(statsList) ):
                drawCompareContiguityPlot( options, statsList[i], statsList[j] )


if __name__ == "__main__":
    main()
