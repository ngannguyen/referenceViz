#!/usr/bin/env python

"""
Create coverage plots
nknguyen at soe dot ucsc dot edu
Input: coverageStats.xml files
"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

#from numpy import *
from numpy import arange
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedFormatter
from matplotlib.font_manager import FontProperties

class Sample(list):
    def __init__(self, sampleNode):
        self.name = sampleNode.attrib[ 'sampleName' ]
        self.referenceName = sampleNode.attrib[ 'referenceName' ]
        self.otherReferenceName = sampleNode.attrib[ 'otherReferenceName' ]
        #self.baseCoverages = sampleNode.attrib[ 'baseCoverages' ].strip().split()
        items = sampleNode.attrib[ 'baseCoverages' ].strip().split()
        totalBases = 0
        for i in range( len(items) ):
            items[i] = int( items[i] )
            totalBases += items[i]
        self.baseCoverages = items

        self.relativeBaseCoverages = []
        for item in items:
            self.relativeBaseCoverages.append( item/float(totalBases) )
        
        self.referenceBasesMapped = int( sampleNode.attrib['referenceBasesMapped'] )
        self.otherReferenceBasesMapped = int( sampleNode.attrib['otherReferenceBasesMapped'] )
        self.totalBases = totalBases
        #list = []
        #culm = 0
        #for i in range( len(items) - 1, -1, -1 ):
        #    culm += items[i]
        #    list.append( culm/float(totalBases) )
        #self.culmBaseCoverages = list

class Stats( list ): #each Stats represents one input XML file
    def __init__( self, name ):
        self.name = name
    def setRefName( self, refname ):
        self.refname = refname
    def setOtherRefName( self, name ):
        self.otherRefname = name

def drawData( axes, stats, isAbs ):
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = []
    ydataList = [] 
    #initialize ydataList:
    #for i in range( len(stats[0].baseCoverages) - len(stats), len(stats[0].baseCoverages) ):
    for i in range( len(stats) -1 ): #each coverage level
        ydata = []
        for j in range( len(stats) ):#each sample
            if isAbs:
                if stats[j].name == 'aggregate':
                    ydata.append( stats[j].baseCoverages[i]/(len(stats) -1) )
                else:
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
    l = axes.bar( x, ydataList[ len(ydataList) - 1 ], barwidth, color = colors[colorindex], ec="w" ) 
    lines.append( l[0] )
    linenames.append( "%d"  % len(ydataList) )
    culmulativeList = ydataList[ len(ydataList) - 1 ]

    for i in range( len(ydataList) - 2, -1, -1 ):
        colorindex += 1
        l = axes.bar( x, ydataList[i], barwidth, color = colors[colorindex], bottom=culmulativeList, ec="w" )
        lines.append( l[0] )
        linenames.append( "%d" % (i + 1) )
        
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
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    if not isAbs:
        axes.set_ylim(0, 1)
    axes.set_xlim(-0.5, len(stats) )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.8, box.height] )

    lines.reverse()
    linenames.reverse()
    legend = axes.legend( lines, linenames, prop=fontP, loc="best", bbox_to_anchor=(1,0.75) )
    legend._drawFrame=False
    
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
    libplot.writeImage( fig, pdf, options )

#================= Mapped bases of each sample onto the reference vs hg19 ===========
def drawCompareData( axes, options, stats, isAbs ):
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = [ stats.otherRefname, stats.refname, "sampleTotalBases" ]

    barwidth = 0.25
    #X data:
    x1data = arange( len(stats) )
    x2data = x1data + barwidth
    x3data = x2data + barwidth

    if isAbs:
        y1data = [ sample.otherReferenceBasesMapped for sample in stats ]
        y2data = [ sample.referenceBasesMapped for sample in stats ]
        y3data = [ sample.totalBases for sample in stats ]
    else:
        y1data = [ 100.0*sample.referenceBasesMapped/sample.totalBases for sample in stats ]
        y2data = [ 100.0*sample.otherReferenceBasesMapped/sample.totalBases for sample in stats ]
        y3data = [ 100.0*sample.totalBases/sample.totalBases for sample in stats ]
    #Average aggregate data:
    numSamples = len(y1data)
    y1data[ numSamples - 1 ] /= float(numSamples - 1)
    y2data[ numSamples - 1 ] /= float(numSamples - 1)
    y3data[ numSamples - 1 ] /= float(numSamples - 1)

    colors =["#1B9E77", "#D95F02", "#7570B3"] 
    #colors =["#EDF8B1", "#7FCDBB", "#2C7FB8"]
    #colors =["#A1DAB4", "#41B6C4", "#225EA8"]
    l1 = axes.bar( x1data, y1data, barwidth, color = colors[0], ec="w" ) 
    lines.append( l1[0] )
    l2 = axes.bar( x2data, y2data, barwidth, color = colors[1], ec="w" ) 
    lines.append( l2[0] )
    l3 = axes.bar( x3data, y3data, barwidth, color = colors[2], ec="w" )
    lines.append( l3[0] )

    libplot.editSpine( axes )
    axes.set_title("Coverage across samples")
    pyplot.xlabel("Samples")
    pyplot.ylabel("Mapped Bases")

    #set ticks:
    samples = []
    for sample in stats:
        samples.append( sample.name )
    fontP = FontProperties()
    fontP.set_size('small')
    #pyplot.xticks( x + barwidth/2., samples, rotation=45, fontproperties=fontP )
    pyplot.xticks( x2data, samples, rotation=45, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )
    
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    axes.set_ylim( 0, max([max(y1data), max(y2data), max(y3data)]) )
    axes.set_xlim(-0.5, len(stats) )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.95, box.height*0.9] )

    legend = axes.legend( lines, linenames, prop=fontP, loc="best", bbox_to_anchor=(0.2,1) )
    #legend = axes.legend( lines, linenames, prop=fontP, loc="upper left", bbox_to_anchor=(1,0.75) )
    legend._drawFrame=False
    
    return 


def drawCompareCoveragePlot( options, stats, isAbs ):
    prefix = "cmpCoverage_"
    if not isAbs:
        prefix = "cmpRelCoverage_"
    options.out = os.path.join( options.outdir, "%s%s_%s" %(prefix, stats.refname, stats.otherRefname) )
    fig, pdf = libplot.initImage( 12.0, 8.0, options )
    axes = fig.add_axes( [0.09, 0.2, 0.9, 0.6] )

    drawCompareData( axes, options, stats, isAbs )
    libplot.writeImage( fig, pdf, options )

#=================
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
            stats.setOtherRefName( stats[0].otherReferenceName )
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
    
    #HACK to temporarily compensate for ben's coverage bug:
    #for stats in statsList:
    #    prevRefBases = 0
    #    prevOtherRefBases = 0
    #    for i in range( len(stats) -1 ):
    #        stats[i].referenceBasesMapped -= prevRefBases
    #        stats[i].otherReferenceBasesMapped -= prevOtherRefBases
    #        prevRefBases += stats[i].referenceBasesMapped
    #        prevOtherRefBases += stats[i].otherReferenceBasesMapped
    #END HACK

    #Sort the statslist in the increasing order of the singleton coverage:
    for i in range( len(statsList) ):
        statsList[i] = sorted( statsList[i], key=lambda x: x.baseCoverages[0] )
    
    for stats in statsList:
        #for sample in stats:
            #print sample.culmBaseCoverages
        drawCoveragePlot( options, stats, True )
        drawCoveragePlot( options, stats, False )

        #Sort stats by referenceBasesMapped:
        sortedStats = Stats( stats[0].name )
        sortedStats.setRefName( stats[0].referenceName )
        sortedStats.setOtherRefName( stats[0].otherReferenceName )
        sortedStats.extend( sorted( stats, key=lambda s:s.referenceBasesMapped ) )
        for i in range( len(sortedStats) ):
            s = sortedStats[i]
            if s.name == s.referenceName or s.name == s.otherReferenceName:
                break
        if i < len(sortedStats):
            sortedStats.pop(i)
                
        drawCompareCoveragePlot( options, sortedStats, True )

    #if len(statsList) >= 2:
    #    for i in range( len(statsList) -1 ):
    #        for j in range( i+1, len(statsList) ):
    #            drawCompareCoveragePlot( options, statsList[i], statsList[j])

if __name__ == "__main__":
    main()





