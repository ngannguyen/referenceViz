#!/usr/bin/env python

"""
Create coverage plots
nknguyen at soe dot ucsc dot edu
Input: coverageStats.xml files
"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedFormatter
from matplotlib.font_manager import FontProperties

class Sample(list):
    def __init__(self, sampleNode):
        self.name = sampleNode.attrib[ 'sampleName' ]
        self.referenceName = sampleNode.attrib[ 'referenceName' ]
        #self.baseCoverages = sampleNode.attrib[ 'baseCoverages' ].strip().split()
        items = sampleNode.attrib[ 'baseCoverages' ].strip().split()
        sum = 0
        for i in range( len(items) ):
            items[i] = int( items[i] )
            sum += items[i]
        for i in range( 1, len(items) ):
            items[i] += items[i -1] 
        for i in range( len(items) ):
            items[i] /= float(sum)
        self.baseCoverages = items

class Stats( list ): #each Stats represents one input XML file
    def __init__( self, name ):
        self.name = name
    def setRefName( self, refname ):
        self.refname = refname

def drawData( axes, stats ):
    lines = []
    linenames = []
    ydataList = [] 
    #initialize ydataList:
    for i in range( len(stats[0].baseCoverages) ):
        ydata = []
        for j in range( len(stats) ):#each sample
            ydata.append( stats[j].baseCoverages[i] )
        ydataList.append(ydata)
    
    colors = libplot.getColors2( len(stats) )
    colorindex = -1
    for i in range( len(stats) -1, -1, -1 ):
        colorindex += 1
        #l = axes.plot( ydataList[i], color=colors[colorindex], marker='.', markersize=9.0 )
        l = axes.fill_between( x=range(len(ydataList[i])), y1=ydataList[i], y2=[0] * len(ydataList[i]) , facecolor=colors[colorindex], linewidth = 0.0)
        lines.append( l )
        linenames.append( i )
    libplot.editSpine( axes )
    axes.set_title("Coverage across samples")
    pyplot.xlabel("Samples")
    pyplot.ylabel("Cumulative coverage")
    return lines, linenames

def drawCoveragePlot( options, stats ):
    options.out = os.path.join(options.outdir, "coverage_" +  stats.refname)
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = libplot.setAxes( fig )

    lines, linenames = drawData( axes, stats )
    #legend = pyplot.legend( lines, linenames, numpoints=1 )
    #legend._drawFrame=False

    #libplot.setTicks( axes )
    samples = []
    for sample in stats:
        samples.append(sample.name)
    axes.set_xticks( range(0, len(stats)) )
    axes.set_xticklabels( samples )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation( 90 )

    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    libplot.writeImage( fig, pdf, options )


def readfiles( options ):
    statsList = [] #each element represents one input xml file
    for f in options.files:
        name = os.path.basename( f ).split( '.' )[0]
        #print "file %s" %name   
        stats = Stats(name)

        xmltree = ET.parse( f )
        root = xmltree.getroot()
        for sample in root.findall( 'statsForSample' ):
            name = sample.attrib[ 'sampleName' ]
            if name != '' and name != 'ROOT':
                stats.append( Sample( sample ) )
        if len(stats) > 0:
            stats.setRefName( stats[0].referenceName )
        statsList.append( stats )

    return statsList

def initOptions( parser ):
    parser.add_option('--title', dest='title', default='Coverage statistics', help='Based title of the plots, default=%default')
    parser.add_option('--ycutoff', dest='ycutoff', default=0.5, type='float')
    parser.add_option('--outdir', dest='outdir', default='.')

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error('Please specify at least one coverageStat.xml file\n')
    options.files = []
    for f in args:
        if not os.path.exists( f ):
            parser.error( '%s does not exist\n' %f )
        else:
            options.files.append( os.path.abspath( f ) )

def main():
    usage = ( 'usage: %prog [options] file1.xml file2.xml\n\n'
              '%prog takes in coverageStats.xml files and create an image file' )
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()

    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    statsList = readfiles( options )

    #Sort the statslist:
    #statsList = sorted( statsList, key=lambda x: x.baseCoverages[0], reverse=True)

    for stats in statsList:
        drawCoveragePlot( options, stats )

    if len(statsList) >= 2:
        for i in range( len(statsList) -1 ):
            for j in range( i+1, len(statsList) ):
                drawCompareCoveragePlot( options, statsList[i], statsList[j])

if __name__ == "__main__":
    main()





