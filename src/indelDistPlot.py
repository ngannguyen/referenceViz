#!/usr/bin/env python

"""
Create indel distribution plots
nknguyen at soe dot ucsc dot edu
June 16 2011
Input: pathStats_*.xml files
Output: Plots of insertion and deletion distribution
"""
import os, sys, math
from optparse import OptionParser
import xml.etree.ElementTree as ET

from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.font_manager import FontProperties

def getSample(samples, name):
    for s in samples:
        if s.attrib[ 'sampleName' ] == name:
            return s
    return None

def getLeaves( allsamples ):
    #Return only samples that are leaves in the tree
    samples = []
    for sample in allsamples:
        name = sample.attrib[ 'sampleName' ]
        if name != "ROOT" and name != "":
        #if name == 'apd':
            samples.append(sample)
    return samples

def getSampleNames( samples ):
    names = []
    for sample in samples:
        names.append( sample.attrib[ 'sampleName' ] )
    return names

def setAxes( fig, numSamples, samplesPerPlot ):
    numaxes = numSamples / samplesPerPlot 
    if numSamples % samplesPerPlot != 0:
        numaxes += 1
    numaxes *= 2 #display plots of insertion & deletion in the same fig

    axesList = []
    axleft = 0.08
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.08
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.1

    h = ( axheight - ( margin * (numaxes - 1) ) )/float( numaxes )

    bottom = axtop - h
    for i in range( numaxes ):#add axes from top down
        axesList.append( fig.add_axes( [axleft, bottom, axwidth, h] ) )
        bottom = bottom - (margin + h)
    
    return axesList

def getFreq( dist, xlogscale, ylogscale ):
    freqDict = {}
    for v in dist:
        if v in freqDict.keys():
            freqDict[ v ] += 1
        else:
            freqDict[ v ] = 1
    
    xdata = sorted( freqDict )
    ydata = []
    for x in xdata:
        ydata.append( freqDict[ x ] )
    
    if xlogscale == "true":
        xdata = log( array(xdata) )
    if ylogscale == "true":
        ydata = log( array(ydata) )

    return xdata, ydata

def drawData( axesList, samples, samplesPerPlot, options ):
    #colors = libplot.getColors2( len(keys) )
    #markers = [".", "s", "^", "--"]
    if len(axesList) %2 != 0:
        sys.stderr.write( 'Number of axes must be even. Got %d\n' %len(axesList) )
        sys.exit( 1 )

    colors = libplot.getColors0()
    #styles = []

    c = -1
    textsize = 'x-small'
    linesDict = {}
    labelsDict = {}
    for i in range( len(axesList)/2 ):
        inslines = []
        dellines = []
        sampleNames = []
        insAxes = axesList[ i ]
        delAxes = axesList[ i  + len(axesList)/2 ]
        
        startIndex = i * samplesPerPlot
        endIndex = min( [startIndex + samplesPerPlot, len(samples)] )
        for j in range( startIndex, endIndex ):
            sample = samples[j]
            sampleNames.append( sample.attrib[ 'sampleName' ] )
            insDist = [int(val) for val in sample.attrib[ 'insertionSizeDistribution' ].split()]
            insXdata, insYdata = getFreq( insDist, options.xlogscale, options.ylogscale )
            delDist = [int(val) for val in sample.attrib[ 'deletionSizeDistribution' ].split()]
            delXdata, delYdata = getFreq( delDist, options.xlogscale, options.ylogscale )
            
            c += 1
            il = insAxes.plot( insXdata, insYdata, color=colors[c] )
            dl = delAxes.plot( delXdata, delYdata, color=colors[c] )

            inslines.append( il )
            dellines.append( dl )
        
        linesDict[ i ] = inslines
        labelsDict[ i ] = sampleNames
        linesDict[ i + len(axesList)/2 ] = dellines
        labelsDict[ i + len(axesList)/2 ] = sampleNames
        #fontp = FontProperties()
        #fontp.set_size( 'x-small' )
        insAxes.set_title( 'Insertions%d' %i )    
        delAxes.set_title( 'Deletions%d' %i )
        
    for i in range( len(axesList) ):
        axes = axesList[ i ]
        #if options.xlogscale == "true":
        #    axes.set_xscale('log', basex=2)
        #if options.ylogscale == "true":
        #    axes.set_yscale('log', basey=2)
        libplot.editSpine( axes )
        
        
        if options.xlogscale == "true":
            axes.set_xlabel('Log 2 of length (bp)', size = textsize)
        else:
            axes.set_xlabel('Length (bp)', size = textsize)
        
        if options.ylogscale == "true":
            axes.set_ylabel('Log 2 of count', size = textsize)
        else:
            axes.set_ylabel('Count', size = textsize)

        #Legend
        legend = axes.legend( linesDict[ i ], labelsDict[ i ], 'upper right', ncol=3 )
        for t in legend.get_texts():
            t.set_fontsize('x-small')
        legend._drawFrame = False

        for label in axes.get_xticklabels():
            #label.set_rotation(45)
            label.set_fontsize( textsize )
        for label in axes.get_yticklabels():
            label.set_fontsize( textsize )
        #box = axes.get_position()
        #axes.set_position( [box.x0, box.y0, box.width*0.8, box.height] )
        #legend = pyplot.legend( lines, options.keys, numpoints=1, prop=fontP, loc="best", bbox_to_anchor=(1, 0.9) )
        #legend._drawFrame = False

        #libplot.setTicks( axes )
        #axes.set_xticks( range( 0, len(samples) ) )
        #axes.set_xticklabels( sampleNames )
        #for label in axes.xaxis.get_ticklabels():
        #    label.set_rotation( 90 )

        #axes.xaxis.set_ticks_position( 'bottom' )
        #axes.yaxis.set_ticks_position( 'left' )

        #axes.set_xlim( -0.5, len(samples) - 0.5 )
        #axes.set_ylim( -20, 6000 )

    return

def drawPlots( options, samples ):
    #sort samples:
    samples = sorted( samples, key=lambda s: s.attrib[ 'sampleName' ] )

    if len(samples) < 1:
        return

    refname = samples[0].attrib[ 'referenceName' ]
    options.out = os.path.join( options.outdir, 'indelDist_' + refname )
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    
    samplesPerPlot = 8
    axesList = setAxes( fig, len(samples), samplesPerPlot )

    lines = drawData( axesList, samples, samplesPerPlot, options )

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
    parser.add_option('--outdir', dest='outdir', default='.', help='Output directory')
    parser.add_option('--xlog', dest='xlogscale', default="true")
    parser.add_option('--ylog', dest='ylogscale', default="true")

def checkOptions( args, options, parser ):
    if len(args) < 1:
        parser.error('Please specify at least one pathStats file.\n')
    options.files = []
    for file in args:
        if not os.path.exists( file ):
            parser.error('File %s does not exit\n' % file)
        options.files.append( os.path.abspath( file ) )
    #options.keys = options.keys.split(',')

def main():
    usage = ('usage: %prog [options] file1.xml file2.xml\n\n'
             '%prog reads in pathStats.xml files and create insertion and deltion distribution plots')

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
        drawPlots( options, samples )

    #if len( statsListF ) >= 2:
    #    for i in range( len(statsListF) -1 ):
    #        for j in range( i + 1, len(statsListF) ):
    #            drawCompareN50Plot( options, statsListF[i], statsListF[j] )


if __name__ == '__main__':
    main()

