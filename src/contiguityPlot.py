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

#from numpy import *
from numpy import linspace
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

class Bucket:
    def __init__( self, bucketElement ):
        self.start = int( bucketElement.attrib[ 'from' ] )
        self.end = int( bucketElement.attrib[ 'to' ] )
        self.mid = ( self.start + self.end )/2
        self.correct = int( bucketElement.attrib[ 'correct' ] )
        self.samples = int( bucketElement.attrib[ 'samples' ] )
        self.aligned = int( bucketElement.attrib[ 'aligned' ] )
        self.correctPerSample = float( bucketElement.attrib[ 'correctPerSample' ] )
        self.correctPerAligned = float( bucketElement.attrib[ 'correctPerAligned' ] )
        self.cumulativeCorrect = int( bucketElement.attrib[ 'cumulativeCorrect' ] )
        self.cumulativeSamples = int( bucketElement.attrib[ 'cumulativeSamples' ] )
        self.cumulativeAligned = int( bucketElement.attrib[ 'cumulativeAligned' ] )
        self.cumulativeCorrectPerSample = float( bucketElement.attrib[ 'cumulativeCorrectPerSample' ] )
        self.cumulativeCorrectPerAligned = float( bucketElement.attrib[ 'cumulativeCorrectPerAligned' ] )
        
class Sample( list ):
    def __init__(self, name, reference):
        self.name = name
        self.reference = reference

    def setBuckets( self, sampleElement ):
        for bucket in sampleElement.findall( 'bucket' ):
            self.append( Bucket( bucket ) )
    
    def setBuckets2( self, buckets ):
        for bucket in buckets:
            self.append( bucket )

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

def setAxisLimits( axes, ycutoff ):
    axes.set_xscale('log')
    axes.set_ylim( ycutoff, 1.001 )

def drawLegend( axes, lines, sampleNames, options ):
    fontP = FontProperties()
    fontP.set_size('small')
    box= axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    #legend = pyplot.legend( lines, sampleNames, numpoints = 1, prop= fontP, loc="best", bbox_to_anchor=(1, 0.5))
    if not options.legendElements:
        legend = pyplot.legend( lines, sampleNames, prop= fontP, loc="best", bbox_to_anchor=(1,0.5))
        legend._drawFrame=False
    elif len(lines) == len(options.legendElements):
        legend = pyplot.legend( lines, options.legendElements, prop= fontP, loc="best", bbox_to_anchor=(1,0.5) )
        legend._drawFrame=False
    else:
        sys.stderr.write('Number of items in --legendElements is different '
                         'from the number of lines plotted\n' )

def drawData( axes, stats, options ):
    #halfsize = len(stats)/2 + len(stats)%2
    #colors = libplot.getColors2( halfsize )
    #colors = libplot.getColors2( len(stats) )
    #styles = { 0:'-', 1:'--' }

    colors = libplot.getColors0()

    #===========

    #dash = 0
    colorindex = -1
    lines = []
    sampleNames = []

    for sample in stats:
        sampleNames.append(sample.name)
        xdata = []
        ydata = []
        for bucket in sample:
            xdata.append( bucket.mid )
            if options.includeCov:
                ydata.append( bucket.correctPerSample )
            else:
                ydata.append( bucket.correctPerAligned )
        
        #if not dash:
        #    colorindex += 1
        colorindex +=1

        l = axes.plot( xdata, ydata, color=colors[colorindex], linewidth=0.7 )
        #l = axes.plot( xdata, ydata, color=colors[colorindex], linestyle=styles[dash], linewidth=0.5 )
        lines.append(l)
        
        #dash = not dash
    
    libplot.editSpine( axes )
    axes.set_title(options.title)
    pyplot.xlabel("Distance")
    pyplot.ylabel("Correct proportion")
    return lines, sampleNames

def drawAggData( axes, data, sortedIndex, xmin, xmax, cutoff, nbins=10 ):
    data = sorted( data, key=lambda point:point[sortedIndex] )
    #Above the x axis is point[sortedIndex] < point[1-sortedIndex]
    updata = []
    equaldata = []
    downdata = []
    for p in data:
        if p[ 0 ] < cutoff or p[ 1 ] < cutoff:
            continue

        if p[ sortedIndex ] < p[ 1 - sortedIndex ]:
            updata.append( p[ sortedIndex ] )
        elif p[ sortedIndex ] == p[ 1 - sortedIndex ]:
            equaldata.append( p[ sortedIndex ] )
        else:
            downdata.append( p[ sortedIndex ] )

    if sortedIndex == 0:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'
    
    #bins = linspace( xmin, xmax, nbins )
    bins = linspace( cutoff, xmax, nbins )
    ymin, ymax = libplot.bihist( updata, downdata, axes, bins, orientation, color='#0198E1' )
    n3, bins3, patch3 = axes.hist( equaldata, bins=bins, orientation = orientation, color = '#800000' )
    
    if sortedIndex == 0:
        ymax3 = max( [ i.get_height() for i in patch3 ] )
    else:
        ymax3 = max( [ i.get_width() for i in patch3 ] )
    ymax = max( ymax3, ymax )
    return ymin, ymax

def intersect(sample1, sample2):
    #sample1 and sample2 must be sorted by the 'mid' field
    s1 = Sample( sample1.name, sample1.reference )
    s2 = Sample( sample2.name, sample2.reference )

    buckets1 = []
    buckets2 = []
    
    i1 = 0
    i2 = 0
    while i1 < len(sample1) and i2 < len(sample2):
        b1 = sample1[i1]
        b2 = sample2[i2]
        if b1.mid == b2.mid: #sample bucket
            buckets1.append(b1)
            buckets2.append(b2)
            i1 += 1
            i2 += 1
        elif b1.mid > b2.mid:
            i2 += 1
        else:
            i1 += 1

    s1.setBuckets2( buckets1 )
    s2.setBuckets2( buckets2 )

    return s1, s2

def drawCompareData( axesList, xstats, ystats, options ):
    #Only draw the overlapped samples:
    #colors = libplot.getColors2( len(xstats) )
    colors = libplot.getColors0()
    colorindex = -1
    lines = []
    sampleNames = []
    p0axes = axesList[0] #plot 0 axes (see def 'setCompareAxes')
    aggData = [] #data points (buckets) of all samples
    
    for xsample in xstats:
        ysample = getSample( ystats, xsample.name )
        if ysample is None:
            continue
        xsample, ysample = intersect(xsample, ysample)
        #if len(xsample) != len(ysample): 
        #    xsample, ysample = intersect(xsample, ysample)
        #    sys.stderr.write( "Error: Two xml files do not have the same number of buckets for sample %s\n" % xsample.name )
            #sys.exit( 1 )
        
        data = [] #list of (x,y) tuples
        colorindex += 1
        for i in range( len( xsample ) ): #each bucket
            if xsample[i].mid != ysample[i].mid:
                sys.stderr.write( "Two xml files have different buckets\n " )
                sys.exit( 1 )
            if options.includeCov:
                data.append( (xsample[i].correctPerSample, ysample[i].correctPerSample) )
            else:
                data.append( (xsample[i].correctPerAligned, ysample[i].correctPerAligned) )

        x2data = [ point[0] for point in data ]
        y2data = [ point[1] for point in data ]
        l = p0axes.plot( x2data, y2data, color=colors[colorindex], marker='.', markersize=4.0, linestyle='none' )
        lines.append( l )
        sampleNames.append( xsample.name )
        aggData.extend( data )

    #Draw the y=x line
    x = [0, 1]
    y = [0, 1]
    p0axes.plot(x, y, color="#919191")

    libplot.editSpine( p0axes )
    p0axes.set_title(options.title)
    p0axes.set_xlabel(xstats.refname)
    p0axes.set_ylabel(ystats.refname)
    libplot.setTicks( p0axes )
    
    #legend:
    fontP = FontProperties()
    fontP.set_size('xx-small')
    legend = p0axes.legend( lines, sampleNames, 'lower right', numpoints = 1, prop=fontP, ncol = 2)
    legend._drawFrame = False
    
    #p0axes.set_xlim( -0.005, 1.005 )
    #p0axes.set_ylim( -0.005, 1.005 )
    p0axes.set_xlim( options.ycutoff, 1 + (1 - options.ycutoff)*0.01 )
    p0axes.set_ylim( options.ycutoff, 1 + (1 - options.ycutoff)*0.01 )
   
    #box = p0axes.get_position()
    #p0axes.set_position([box.x0, box.y0, box.width * 0.8, box.height * 0.8])
    #legend = pyplot.legend( lines, sampleNames, numpoints = 1, prop= fontP, loc="best", bbox_to_anchor=(1, 0.6))
    #legend._drawFrame=False
    
    #DRAW AGGREGATE DATA (plot 1 and plot 2):
    nbins = 20
    p1axes = axesList[1]
    y1min, y1max = drawAggData( p1axes, aggData, 0, 0, 1, options.ycutoff, nbins )
    y1lim = max( abs(y1min), abs(y1max) )
    p1axes.set_ylim( -y1lim*1.1, y1lim*1.1 )
    p1axes.set_xlim( options.ycutoff, 1 + (1-options.ycutoff)*0.01 )
    #p1axes.set_ylim( y1min*1.1, y1max*1.1 )
    for loc, spine in p1axes.spines.iteritems():
        if loc == 'left':
            spine.set_position( ( 'outward', 10 ) )
        spine.set_color( 'none' )
    p1axes.axhline( 0, color = '#000000' )
    p1axes.xaxis.set_major_locator( NullLocator() )
    p1axes.xaxis.set_major_formatter( NullFormatter() )
    p1axes.yaxis.set_ticks([-y1lim, 0, y1lim])
    #p1axes.yaxis.set_ticklabels([-y1lim, 0, y1lim])
    for l in p1axes.yaxis.get_ticklabels():
        l.fontproperties = fontP
    #p1axes.tick_params( axis='y', labelsize='small')

    p2axes = axesList[2]
    x2min, x2max = drawAggData( p2axes, aggData, 1, 0, 1, options.ycutoff, nbins )
    x2lim = max( abs(x2min), abs(x2max) )
    p2axes.set_xlim( -x2lim*1.1, x2lim*1.1 )
    p2axes.set_ylim( options.ycutoff, 1 + (1- options.ycutoff)*0.01 )
    #p2axes.set_xlim( x2min*1.1, x2max*1.1 )
    for loc, spine in p2axes.spines.iteritems():
        if loc == 'bottom':
            spine.set_position( ( 'outward', 10 ) )
        spine.set_color( 'none' )
    p2axes.axvline( 0, color = '#000000' )
    p2axes.yaxis.set_major_locator( NullLocator() )
    p2axes.yaxis.set_major_formatter( NullFormatter() )
    p2axes.xaxis.set_ticks([-x2lim, 0, x2lim])
    for l in p2axes.xaxis.get_ticklabels():
        l.fontproperties = fontP
        l.set_rotation( 45 )
    return


def drawContiguityPlot( options, stats ):
    #options.out = os.path.join(options.outdir, "contiguity_" + stats.refname) #name of output file
    options.out = os.path.join(options.outdir, options.exp + "_" + stats.refname) #name of output file
    if options.includeCov:
        options.out = options.out + "_incCov"
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = libplot.setAxes( fig )
    
    lines, sampleNames = drawData( axes, stats, options )
    drawLegend( axes, lines, sampleNames, options )
    setAxisLimits( axes, options.ycutoff )
    libplot.setTicks( axes )

    libplot.writeImage( fig, pdf, options )

def setCompareAxes( fig ):
    """
    Set axes for the CompareContiguityPlot. There are 3 subplots total:
    Plot 0: The main plot where each axis represents one reference seq. Eg: xaxis <- cactusref, yaxis <- hg19
    Plot 1: Right at the top of plot 1, shows frequencies of data points that lie above and below the y=x axis
    Plot 2: To the right of plot 1, shows frequencies of data points that lie on the left and right the y=x axis
    """
    axesList = []
    axleft = 0.1
    axright = 0.98
    axwidth = axright - axleft
    axbottom = 0.08
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.07 #space between plots

    plot0wh = 0.7 #plot1 width = plot1 height
    plot1h = axheight - (plot0wh + margin) 
    plot2w = plot1h

    axesList.append( fig.add_axes( [axleft, axbottom, plot0wh, plot0wh] ) ) #Plot 0
    axesList.append( fig.add_axes( [axleft, axbottom + plot0wh + margin, plot0wh, plot1h] ) ) #Plot 1
    axesList.append( fig.add_axes( [axleft + plot0wh + margin, axbottom, plot2w, plot0wh] ) ) #Plot 2
    return axesList

def drawCompareContiguityPlot( options, xstats, ystats ):
    #options.out = os.path.join(options.outdir, "contiguity_" + xstats.refname + "_" + ystats.refname)
    options.out = os.path.join(options.outdir, options.exp + "_" + xstats.refname + "_" + ystats.refname)
    if options.includeCov:
        options.out = options.out + "_incCov"
    fig, pdf = libplot.initImage( 8.0, 8.0, options )
    
    #Set axes:
    #axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )
    axesList = setCompareAxes( fig )

    drawCompareData( axesList, xstats, ystats, options )
    
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
                s = Sample( name, sample.attrib[ 'referenceName' ] )
                s.setBuckets( sample )
                stats.append( s )
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
    parser.add_option('--outdir', dest='outdir', default='.', help='Output directory')
    parser.add_option('--includeCoverage', dest='includeCov', action="store_true", default=False, help='If specified, will include coverage info in the plots')

def checkOptions( args, options, parser ):
    options.files = []
    for f in args:
        if not os.path.exists( f ):
            parser.error('%s does not exist\n' %f)
        options.files.append( os.path.abspath( f ) )
    
    if len(options.files) < 1:
        parser.error('Please specify at least one valid contiguityStats file.\n')
    options.exp = ( os.path.basename( options.files[0] ).split('_') )[0]

    #system("mkdir -p %s" % options.outdir)
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
