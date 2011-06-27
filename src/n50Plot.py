#!/usr/bin/env python

"""
Create N50 plots
nknguyen at soe dot ucsc dot edu
June 14 2011
Input:

"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.font_manager import FontProperties

#class Sample():
#    def __init__( self, sampleElement ):
#        self.name = sampleElement.attrib[ 'sampleName' ]
#        self.refname = sampleElement.attrib[ 'referenceName' ]
#        self.
def drawN50data( axes, samples, options ):
    #key can be 'blockN50', 'sequenceN50', 'contigPathN50', or 'scaffolfPathN50'
    keys = options.keys
    colors = libplot.getColors0()
    #markers = [".", "s", "^", "--"]

    c = -1
    lines = []
    for key in keys:
        ydata = []
        for sample in samples:
            y = int(sample.attrib[ key ])
            ydata.append( y )
        if options.logscale:
            ydata = log10( array(ydata) )
        
        c += 1
        l = axes.plot( ydata, color=colors[c], marker=".", markersize=10.0, linestyle='none' )
        lines.append(l)

    libplot.editSpine( axes )
    pyplot.xlabel( 'Samples' )
    if options.logscale:
        pyplot.ylabel( 'Log 10 of N50' )
    else:
        pyplot.ylabel( 'N50' )

    return lines

def getSample(samples, name):
    for s in samples:
        if s.attrib[ 'sampleName' ] == name:
            return s
    return None

def drawCompareN50data( axes, xsamples, ysamples, options ):
    keys = options.keys
    colors = libplot.getColors0()
    c = -1
    lines = []
    xrefname = xsamples[0].attrib[ 'referenceName' ]
    yrefname = ysamples[0].attrib[ 'referenceName' ]

    minval = inf
    maxval = 0
    for key in keys:
        xdata = []
        ydata = []
        for xsample in xsamples:
            name = xsample.attrib[ 'sampleName' ]
            if name == yrefname:
                continue
            ysample = getSample( ysamples, name )
            if ysample == None:
                sys.stderr.write( "%s has %s sample, but %s doesn't\n" % (xrefname, name, yrefname) )
                continue
            
            xval = int(xsample.attrib[key])
            yval = int(ysample.attrib[key])
            if xval > 0 and yval > 0:
                xdata.append(xval)
                ydata.append(yval)
            #xdata.append( int(xsample.attrib[ key ]) )
            #ydata.append( int(ysample.attrib[ key ]) )
        
        if options.logscale:
            xdata = log10( array(xdata) )
            ydata = log10( array(ydata) )
         
        c += 1
        l = axes.plot( xdata, ydata, color=colors[c], marker=".", markersize=10.0, linestyle='none' )
        lines.append(l)
        
        currmax = max( xdata.max(), ydata.max() )
        if maxval < currmax:
            maxval = currmax
        
        currmin = min( xdata.min(), ydata.min() )
        if minval > currmin:
            minval = currmin
    
    if minval == -inf:
        minval = 0
    #Draw y=x line
    span = maxval - minval
    #print "MaxVal: %f, MinVal: %f. Span: %f" % (maxval, minval, span)
    x = [ minval - span*0.1, maxval + span*0.1 ]
    y = [ minval - span*0.1, maxval + span*0.1 ]
    axes.plot( x, y, color="0.9" )

    libplot.editSpine( axes )
    pyplot.xlabel( xrefname )
    pyplot.ylabel( yrefname )

    return lines, maxval, minval

def getLeaves( allsamples ):
    #Return only samples that are leaves in the tree
    samples = []
    for sample in allsamples:
        name = sample.attrib[ 'sampleName' ]
        if name != "ROOT" and name != "":
            samples.append(sample)
    return samples

def getSampleNames( samples ):
    names = []
    for sample in samples:
        names.append( sample.attrib[ 'sampleName' ] )
    return names

def drawN50Plot( options, samples ):
    #sort samples:
    samples = sorted( samples, key=lambda s: int(s.attrib[ options.sortkey ]), reverse=True )
    sampleNames = getSampleNames( samples )

    if len(samples) < 1:
        return

    refname = samples[0].attrib[ 'referenceName' ]
    options.out = os.path.join( options.outdir, 'n50_' + refname )
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    title = "N50"
    lines = drawN50data( axes, samples, options )
    axes.set_title(title)

    #Legend
    fontP = FontProperties()
    fontP.set_size( 'small' )
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.8, box.height] )
    legend = pyplot.legend( lines, options.keys, numpoints=1, prop=fontP, loc="best", bbox_to_anchor=(1, 0.9) )
    legend._drawFrame = False

    #libplot.setTicks( axes )
    axes.set_xticks( range( 0, len(samples) ) )
    axes.set_xticklabels( sampleNames )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation( 90 )

    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    axes.set_xlim( -0.5, len(samples) - 0.5 )
    #axes.set_ylim( -20, 6000 )

    libplot.writeImage( fig, pdf, options )

def drawCompareN50Plot( options, xsamples, ysamples ):
    #Sort xsamples and ysamples in the order of the sampleNames:
    xsamples = sorted( xsamples, key=lambda s: s.attrib[ 'sampleName' ] )
    
    if len(xsamples) < 1 or len(ysamples) < 1:
        return

    xrefname = xsamples[0].attrib[ 'referenceName' ]
    yrefname = ysamples[0].attrib[ 'referenceName' ]
    options.out = os.path.join( options.outdir, 'n50_' + xrefname + '_' + yrefname )
    fig, pdf = libplot.initImage( 8.0, 8.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    lines, maxval, minval = drawCompareN50data( axes, xsamples, ysamples, options )
    title = "N50 %s versus %s" % (xrefname, yrefname)
    axes.set_title(title)
     
    #Legend
    fontP = FontProperties()
    fontP.set_size( 'small' )
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.8, box.height*0.8] )
    legend = pyplot.legend( lines, options.keys, numpoints=1, prop=fontP, loc="best", bbox_to_anchor=(1, 0.5) )
    legend._drawFrame = False

    libplot.setTicks( axes )
    span = maxval - minval
    axes.set_xlim( minval - span*0.1, maxval + span*0.1 )
    axes.set_ylim( minval - span*0.1, maxval + span*0.1 )
    libplot.writeImage( fig, pdf, options )


def readfiles( files ):
    statsList = [] #each element represents each xml file
    for f in files:
        xmltree = ET.parse( f )
        statsList.append( xmltree )

    return statsList

def initOptions( parser ):
    #parser.add_option('--title', dest='title', default='',
    #                  help='')
    #parser.add_option('--ycutoff', dest='ycutoff', default=0.5, type='float',
    #                   help='Only points with y values from ycutoff to 1 are displayed')
    parser.add_option('--sortkey', dest='sortkey', default='scaffoldPathN50',
                       help='key attribute for sorting')
    parser.add_option('--keys', dest='keys', default='blockN50,contigPathN50,scaffoldPathN50',
                       help='comma separted list of which N50 statistics to be displayed')
    parser.add_option('--outdir', dest='outdir', default='.', help='Output directory')
    parser.add_option('--abs', dest='logscale', action="store_false", default=True, help="If specified, the axes have the absolute scale, instead of log 10 as default")

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error('Please specify at least one pathStats file.\n')
    options.files = []
    for file in args:
        if not os.path.exists( file ):
            parser.error('File %s does not exit\n' % file)
        options.files.append( os.path.abspath( file ) )
    options.keys = options.keys.split(',')

def main():
    usage = ('usage: %prog [options] file1.xml file2.xml\n\n'
             '%prog reads in pathStats.xml files and create N50 plots')

    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser)
    libplot.checkOptions( options, parser )

    statsList = readfiles( options.files )

    #Filter out non-leaf samples:
    statsListF = []
    for xmltree in statsList:
        root = xmltree.getroot()
        allsamples = root.findall( 'statsForSample' )

        samples = getLeaves( allsamples )
        statsListF.append( samples )

    for samples in statsListF:
        drawN50Plot( options, samples )

    if len( statsListF ) >= 2:
        for i in range( len(statsListF) -1 ):
            for j in range( i + 1, len(statsListF) ):
                drawCompareN50Plot( options, statsListF[i], statsListF[j] )


if __name__ == '__main__':
    main()

