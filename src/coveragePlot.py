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

class Sample():
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
        if totalBases == 0:
            sys.stderr.write("Total bases is 0 for sample %s. Please check.\n" %self.name)
            sys.exit(1)
        self.baseCoverages = items

        self.relativeBaseCoverages = []
        for item in items:
            self.relativeBaseCoverages.append( item/float(totalBases) )
        
        self.referenceBasesMapped = int( sampleNode.attrib['referenceBasesMapped'] )
        self.otherReferenceBasesMapped = int( sampleNode.attrib['otherReferenceBasesMapped'] )
        self.totalBases = totalBases

        #set order: 
        self.order = 0
        orders = {'hg19':1,'reference':2, 'average':3, 'all':4, 'panTro3':6, 'minusOtherReference':5 }

        if self.name in orders:
            self.order = orders[self.name]

    def __cmp__(self, other):
        if self.order > other.order:
            return 1
        elif self.order < other.order:
            return -1
        else:
            if self.baseCoverages[0] > other.baseCoverages[0]:
                return 1
            elif self.baseCoverages[0] < other.baseCoverages[0]:
                return -1
            else:
                return 0

#class Stats( list ): #each Stats represents one input XML file
#    def __init__( self, name ):
#        self.name = name
#    def setRefName( self, refname ):
#        self.refname = refname
#    def setOtherRefName( self, name ):
#        self.otherRefname = name

def drawData( axes, stats, isAbs, ycutoff ):
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = []
    ydataList = [] 
    #initialize ydataList:
    #for i in range( len(stats[0].baseCoverages) - len(stats), len(stats[0].baseCoverages) ):
    #for i in range( len(stats) -1 ): #each coverage level
    for i in range( len(stats) -1 - 2 ): #each coverage level (num samples - average, reference, minusOtherReference
        ydata = []
        for j in range( len(stats) ):#each sample
            if isAbs:
                #if stats[j].name == 'aggregate':
                #    ydata.append( stats[j].baseCoverages[i]/(len(stats) -1) )
                #else:
                ydata.append( stats[j].baseCoverages[i] )
            else:
                ydata.append( stats[j].relativeBaseCoverages[i] )

        ydataList.append(ydata)

    #colors = libplot.getColors2( len(stats) )
    colors = libplot.getColors3()
    colorindex = 0
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
    #axes.set_title("") #TO BE NAMED!!!
    pyplot.xlabel("Samples")
    pyplot.ylabel("Proportion of total bases")

    #set ticks:
    samples = []
    for sample in stats:
        samples.append( libplot.properName( sample.name ) )
    fontP = FontProperties()
    fontP.set_size('small')
    pyplot.xticks( x + barwidth/2., samples, rotation=90, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )    
    
    #for label in axes.yaxis.get_ticklabels():
    #    label.fontproperties = fontP
    #    label.set_rotation( 45 )
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    miny = ycutoff
    if not isAbs:
        axes.set_ylim(ycutoff, 1)
        #axes.set_ylim(0, 1)
    axes.set_xlim(-0.5, len(stats) )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.8, box.height] )

    lines.reverse()
    linenames.reverse()
    legend = axes.legend( lines, [libplot.properName(n) for n in linenames], prop=fontP, loc="best", bbox_to_anchor=(1,0.75) )
    legend._drawFrame=False
    
    return lines, linenames

def drawCoveragePlot( options, stats, isAbs, ycutoff ):
    prefix = "coverage_%.2f_" %ycutoff
    if not isAbs:
        prefix = "rel_coverage_%.2f_" %ycutoff
    options.out = os.path.join(options.outdir, prefix +  stats[0].referenceName)
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.14, 0.2, 0.8, 0.6] )
    #axes = libplot.setAxes( fig )

    lines, linenames = drawData( axes, stats, isAbs, ycutoff )
    libplot.writeImage( fig, pdf, options )

#================= Mapped bases of each sample onto the reference vs hg19 ===========
def drawCompareData( axes, options, stats, isAbs ):
    if len(stats) == 0:
        return
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = [ stats[0].otherReferenceName, stats[0].referenceName, "total" ]

    barwidth = 0.25
    #X data:
    x1data = []
    #avgIndex = -1

    currx = -1
    for i,s in enumerate( stats ):
        #if s.name == 'average':
        #    avgIndex = i
        if s.name == 'average' or s.name == 'panTro3':
            currx += 1 + 1.5*barwidth
        else:
            currx += 1
        x1data.append( currx )

    #print x1data
    x2data = [ x + barwidth for x in x1data ]
    x3data = [ x + barwidth for x in x2data ]

    if isAbs:
        y1data = [ sample.otherReferenceBasesMapped for sample in stats ]
        y2data = [ sample.referenceBasesMapped for sample in stats ]
        y3data = [ sample.totalBases for sample in stats ]
    else:
        y1data = [ 100.0*sample.referenceBasesMapped/sample.totalBases for sample in stats ]
        y2data = [ 100.0*sample.otherReferenceBasesMapped/sample.totalBases for sample in stats ]
        y3data = [ 100.0*sample.totalBases/sample.totalBases for sample in stats ]
    
    #Average aggregate data:
    #if avgIndex > 0:
    #    y1data[ avgIndex ] /= float(avgIndex)
    #    y2data[ avgIndex ] /= float(avgIndex)
    #    y3data[ avgIndex ] /= float(avgIndex)

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
    #axes.set_title("Coverage across samples") #TO BE NAMED
    pyplot.xlabel("Samples")
    pyplot.ylabel("Number of bases")

    #set ticks:
    samples = []
    for sample in stats:
        samples.append( libplot.properName(sample.name) )
    fontP = FontProperties()
    fontP.set_size('small')
    #pyplot.xticks( x + barwidth/2., samples, rotation=45, fontproperties=fontP )
    pyplot.xticks( x2data, samples, rotation=45, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )
    
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    miny = min( [min(y1data), min(y2data), min(y3data)] )
    miny = miny*0.9
    axes.set_ylim( miny, max([max(y1data), max(y2data), max(y3data)]) )
    #axes.set_xlim(-0.5, len(stats) + 0.5 )
    axes.set_xlim(-0.5, max(x3data) + 0.5 )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.95, box.height*0.9] )

    #legend = axes.legend( lines, linenames, prop=fontP, loc="best", bbox_to_anchor=(0.2,1) )
    legend = axes.legend( lines, [libplot.properName(n) for n in linenames], prop=fontP, loc="best", bbox_to_anchor=(0.2, 1) )
    #legend = axes.legend( lines, linenames, prop=fontP, loc="upper left", bbox_to_anchor=(1,0.75) )
    legend._drawFrame=False
    
    return 

def drawCompareCoveragePlot( options, stats, isAbs ):
    if len(stats) == 0:
        return
    prefix = "cmpCoverage_"
    if not isAbs:
        prefix = "cmpRelCoverage_"
    options.out = os.path.join( options.outdir, "%s%s_%s" %(prefix, stats[0].referenceName, stats[0].otherReferenceName) )
    fig, pdf = libplot.initImage( 12.0, 8.0, options )
    axes = fig.add_axes( [0.09, 0.2, 0.9, 0.6] )

    drawCompareData( axes, options, stats, isAbs )
    libplot.writeImage( fig, pdf, options )

#================ display data in scatter plot form, xaxis = number of sample covered, yaxis = number of bases
def drawScatter( axes, options, stats ):
    if len(stats) < 4:
        return
    
    #samples = ["panTro3", "minusOtherReference", "average", "reference", "hg19"]
    samples = ["reference", "minusOtherReference", "hg19", "panTro3", "average"]
    xdata = range( 0, len(stats) -4 )
    #print xdata
    ydataList = []
    miny = float('inf')
    maxy = float('-inf')
    for name in samples:
        for s in stats:
            if s.name == name:
                ydataList.append( s.baseCoverages[: len(stats) -4] )
                miny = min( [miny, min(s.baseCoverages[: len(stats) -4])] )
                maxy = max( [maxy, max(s.baseCoverages[: len(stats) -4])] )
                break

    lines = []
    colors = libplot.getColors6()
    c = -1
    offset = 0.12
    axes.set_yscale('log')
    minorLocator = LogLocator(subs=[1,2,3,4,5])
    axes.yaxis.set_minor_locator(minorLocator)

    for i in xrange( len(samples) ):
        xdatai = [x + offset*i for x in xdata]
        ydata = ydataList[i]
        c += 1
        l = axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=12.0, linestyle='none')
        lines.append(l)
    
    fontP = FontProperties()
    fontP.set_size('x-small')

    yrange = maxy - miny
    miny = miny - yrange*0.15
    maxy = maxy + yrange*0.15
    
    #Draw vertical lines to separate each category:
    for i in xrange(1, len(stats) -4):
        d = (1 - offset*len(samples))/2.0
        x = [i -d, i -d]
        y = [ max([0,miny]), maxy ]
        axes.plot(x, y, color="#CCCCCC", linewidth=0.005)

    xmin = -0.4
    xmax = len(stats) - 4 -1 + offset*len(samples) + offset
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(miny, maxy)
    libplot.editSpine(axes)
    
    axes.set_xticks( [ i + offset*(len(samples)/2-1) for i in range(0, len(stats) -4)] )
    axes.set_xticklabels( range(1, len(stats) -2) )
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    legend = pyplot.legend( lines, [libplot.properName(s) for s in samples], numpoints=1, loc='best', prop=fontP )
    legend._drawFrame = False

    axes.set_xlabel( 'Number of samples' )
    axes.set_ylabel( 'Number of bases' )
    #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    #axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    return

def drawScatterPlot( options, stats ):
    prefix = "coverageScatter_"
    options.out = os.path.join( options.outdir, "%s" %(prefix) )
    fig, pdf = libplot.initImage( 12.0, 8.0, options )
    axes = fig.add_axes( [0.1, 0.15, 0.85, 0.75] )
    drawScatter( axes, options, stats )
    libplot.writeImage( fig, pdf, options )

#=================
def readfiles( options ):
    statsList = [] #each element represents one input xml file
    for f in options.files:
        name = os.path.basename( f ).split( '.' )[0]
        #print "file %s" %name   
        #stats = Stats(name)
        stats = []

        xmltree = ET.parse( f )
        root = xmltree.getroot()
        for sample in root.findall( 'statsForSample' ):
            name = sample.attrib[ 'sampleName' ]
            if name != '' and name != 'ROOT' and name not in options.filteredSamples:
                stats.append( Sample( sample ) )
        #if len(stats) > 0:
        #    stats.setRefName( stats[0].referenceName )
        #    stats.setOtherRefName( stats[0].otherReferenceName )
        statsList.append( stats )

    return statsList

def initOptions( parser ):
    parser.add_option('--title', dest='title', default='Coverage statistics', help='Based title of the plots, default=%default')
    parser.add_option('--ycutoff', dest='ycutoff', default=0.9, type='float')
    parser.add_option('--outdir', dest='outdir', default='.')
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    #parser.add_option('--samplesOrder', dest='samplesOrder', default='', help='Order of the samples to display')

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error('Please specify at least one coverageStat.xml file\n')
    options.files = []
    for f in args:
        if not os.path.exists( f ):
            parser.error( '%s does not exist\n' %f )
        else:
            options.files.append( os.path.abspath( f ) )
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split('-')
    else:
        options.filteredSamples = []

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
    
    for stats in statsList:
        stats.sort()
        stats1 = []
        for s in stats:
            if s.name != 'all':
                stats1.append(s)
        drawCoveragePlot( options, stats1, True, 0 )
        drawCoveragePlot( options, stats1, False, 0 )
        if options.ycutoff > 0:
            drawCoveragePlot( options, stats1, True, options.ycutoff )
            drawCoveragePlot( options, stats1, False, options.ycutoff )

        #sort by totalBases:
        specialcases = {'average':None, 'all':None, 'panTro3':None}
        sortedstats = []
        for i in xrange( len(stats) ):
            if stats[i].name in [ stats[i].referenceName, stats[i].otherReferenceName, 'minusOtherReference' ]:
                continue
            if stats[i].name in specialcases:
                specialcases[ stats[i].name ] = stats[i]
            else:
                sortedstats.append( stats[i] )
        sortedstats = sorted( sortedstats, key=lambda s: s.totalBases )
        for k in specialcases:
            s = specialcases[k]
            if s:
                sortedstats.append( s )
        if len(sortedstats) > 0:
            drawCompareCoveragePlot( options, sortedstats, True )
        
        drawScatterPlot( options, stats )


if __name__ == "__main__":
    main()





