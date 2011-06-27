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
        self.baseCoverages = items

        self.relativeBaseCoverages = []
        for item in items:
            self.relativeBaseCoverages.append( item/float(sum) )

        #list = []
        #culm = 0
        #for i in range( len(items) - 1, -1, -1 ):
        #    culm += items[i]
        #    list.append( culm/float(sum) )
        #self.culmBaseCoverages = list

class Stats( list ): #each Stats represents one input XML file
    def __init__( self, name ):
        self.name = name
    def setRefName( self, refname ):
        self.refname = refname

def drawData( axes, stats, isAbs ):
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = []
    ydataList = [] 
    #initialize ydataList:
    #for i in range( len(stats[0].baseCoverages) - len(stats), len(stats[0].baseCoverages) ):
    for i in range( len(stats) ): #each coverage level
        ydata = []
        for j in range( len(stats) ):#each sample
            if isAbs:
                ydata.append( stats[j].baseCoverages[i] )
            else:
                ydata.append( stats[j].relativeBaseCoverages[i] )

        ydataList.append(ydata)

    #colors = libplot.getColors2( len(stats) )
    colors = libplot.getColors0()
    colorindex = 0
    #for i in range( len(stats) -1, -1, -1 ):
    x = arange( len(stats) ) #x axis represents the samples
    barwidth = 0.6

    #add bottom-most bar (number of bases that are in all samples)
    l = axes.bar( x, ydataList[ len(ydataList) - 1 ], barwidth, color = colors[colorindex] ) 
    lines.append( l )
    culmulativeList = ydataList[ len(ydataList) - 1 ]

    for i in range( len(ydataList) - 2, -1, -1 ):
        colorindex += 1
        l = axes.bar( x, ydataList[i], barwidth, color = colors[colorindex], bottom=culmulativeList )
        lines.append( l )
        linenames.append( i )
        
        #Update cumulative list:
        for j in range( len(culmulativeList) ):
            culmulativeList[j] += ydataList[i][j]
        #l = axes.fill_between( x=range(len(ydataList[i])), y1=ydataList[i], y2=[0] * len(ydataList[i]) , facecolor=colors[colorindex], linewidth = 0.0)
    libplot.editSpine( axes )
    axes.set_title("Coverage across samples")
    pyplot.xlabel("Samples")
    pyplot.ylabel("Coverage")

    #set ticks:
    samples = []
    for sample in stats:
        samples.append( sample.name )
    fontP = FontProperties()
    fontP.set_size('small')
    pyplot.xticks( x + barwidth/2., samples, rotation=45, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )    
    
    #for label in axes.yaxis.get_ticklabels():
    #    label.fontproperties = fontP
    #    label.set_rotation( 45 )
    box = axes.get_position()
    #axes.set_position( [box.x0, box.y0, box.width, box.height * 0.8] )

    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    if not isAbs:
        axes.set_ylim(0, 1)
    axes.set_xlim(-0.5, len(stats) )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    return lines, linenames

def drawCoveragePlot( options, stats, isAbs ):
    prefix = "coverage_"
    if not isAbs:
        prefix = "rel_coverage_"
    options.out = os.path.join(options.outdir, prefix +  stats[0].referenceName)
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.14, 0.2, 0.8, 0.6] )
    #axes = libplot.setAxes( fig )

    lines, linenames = drawData( axes, stats, isAbs )
    #legend = pyplot.legend( lines, linenames, numpoints=1 )
    #legend._drawFrame=False

    #libplot.setTicks( axes )
    #samples = []
    #for sample in stats:
    #    samples.append(sample.name)
    #axes.set_xticks( range(0, len(stats)) )
    #axes.set_xticklabels( samples )
    #for label in axes.xaxis.get_ticklabels():
    #    label.set_rotation( 90 )

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

    #Sort the statslist in the increasing order of the singleton coverage:
    for i in range( len(statsList) ):
        statsList[i] = sorted( statsList[i], key=lambda x: x.baseCoverages[0] )
    
    for stats in statsList:
        #for sample in stats:
            #print sample.culmBaseCoverages
        drawCoveragePlot( options, stats, True )
        drawCoveragePlot( options, stats, False )

    #if len(statsList) >= 2:
    #    for i in range( len(statsList) -1 ):
    #        for j in range( i+1, len(statsList) ):
    #            drawCompareCoveragePlot( options, statsList[i], statsList[j])

if __name__ == "__main__":
    main()





