#!/usr/bin/env python

"""
Create copy number variation plots
nknguyen at soe dot ucsc dot edu
Aug 15 2011
Input: copyNumberStats.xml
Output: cnv plots, one for each sample
"""

import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

import libPlotting as libplot
import matplotlib.pylab  as pylab
import matplotlib.pyplot as pyplot
import matplotlib.patches as patches
import matplotlib.lines as lines
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties


def drawOneCnvPlot( refCn, ax, cnvDict, options, minCn, maxCn ):
    h = 0.2
    l = 0.2
    options.facetText = l
    margin = 0.05

    maxColCount = float( max( cnvDict.values() ) )
    if maxColCount == 0:
        return
        #print refCn
        #print cnvDict
        #sys.exit(1)

    for i in range(minCn, maxCn + 1):
        if str(i) == refCn:
            color = 'k'
            weight = 'medium'
        else:
            color = (0.5, 0.5, 0.5)
            weight = 'normal'
        ax.text( x=l, y=i+0.5, s=str(i), horizontalalignment='right',
                 verticalalignment='center', fontsize=9, color=color, weight=weight)
        if str(i) in cnvDict.keys():
            if str(i) == refCn:
                color = (0.3, 0.3, 0.3)
            else:
                color = 'red'
            ax.add_patch( patches.Rectangle( xy=(l+margin, i+0.5 - h/2.0),
                                             width=cnvDict[str(i)]/maxColCount, height=h,
                                             color=color, edgecolor=None, linewidth=0.0) )

def prettyIntStr( strnum ):
    s = ''
    l = len( strnum )
    numIter = l/3
    if l%3 != 0:
        numIter += 1
    for i in range( numIter ):
        end = l - i*3
        start = max(end - 3, 0)
        if i == 0:
            s = strnum[ start: end ]
        else:
            s = strnum[ start: end ] + ',' + s
    return s

def drawAxisLabels( axDict, cnvDict, options, title, maxCn ):
    i = 0
    for refCn in options.sortedAxesNames:
        i = not i
        axDict[ refCn ].set_title( refCn, fontsize = 9 )
        axDict[ refCn ].text( x = options.facetText,
                              y = - maxCn/( 12.0 + 18.0*i ),
                              s = prettyIntStr( str( max( cnvDict[refCn].values() ) ) ),
                              horizontalalignment = 'left',
                              verticalalignment = 'top',
                              fontsize = 8 )

    transAxes = axDict[ 'bg' ].transAxes
    # grey line across top of plot
    axDict[ 'bg' ].add_line( lines.Line2D( xdata = [ options.axleft, options.axright ],
                                           ydata = [ options.axtop, options.axtop ],
                                           color = (0.8, 0.8, 0.8),
                                           linewidth = 0.1, transform = transAxes ) )
    # plot title
    axDict[ 'bg' ].text( x=0.01, y = 0.98, s = title, horizontalalignment='left', verticalalignment='top', fontsize=10, transform=transAxes )

def getSampleData( sample ): 
    minCn = sys.maxint
    maxCn = -sys.maxint
    cnvDict = {} #key = reference copy number, val = {query copy number: columnCount} 
    for cat in sample.findall('copyNumberCategory'):
        r = cat.attrib['referenceCopyNumber']
        q = cat.attrib['assemblyCopyNumber']
        c = int( cat.attrib['columnCount'] )
        #print r, q, c
        if r not in cnvDict:
            cnvDict[ r ] =  {q: c}
        else:
            cnvDict[ r ][ q ] = c

        currMin =  min( int(r), int(q) )
        currMax =  max( int(r), int(q) )
        if minCn > currMin:
            minCn = currMin
        if maxCn < currMax:
            maxCn = currMax

    #print cnvDict
    return cnvDict, minCn, maxCn

def setAxes( fig, refCategories, options ):
    axDict = {}
    #Background axes:
    axDict[ 'bg' ] = fig.add_axes( [ 0.0, 0.0, 1.0, 1.0 ] )
    axDict[ 'bg' ].xaxis.set_major_locator( pylab.NullLocator() )
    axDict[ 'bg' ].yaxis.set_major_locator( pylab.NullLocator() )
    pyplot.box( on=False )

    #Set axes for the categories plots
    options.axleft = 0.01
    options.axright = 0.99
    options.width = options.axright - options.axleft
    options.axbottom = 0.1
    options.axtop = 0.85
    options.axheight = options.axtop - options.axbottom
    margin = 0.05
    w = ( options.width - margin*(len(refCategories) - 1) )/len(refCategories)
    xpos = options.axleft
    refCategories.sort()
    options.sortedAxesNames = refCategories
    for c in refCategories:
        axDict[ c ] = fig.add_axes( [ xpos, options.axbottom, w, options.axheight ] )
        axDict[ c ].xaxis.set_major_locator( pylab.NullLocator() )
        axDict[ c ].yaxis.set_major_locator( pylab.NullLocator() )
        xpos += w + margin
        pyplot.box( on = False )

    return axDict

def setAxisLimits( axDict, minCn, maxCn ):
    for refCn in axDict:
        axDict[ refCn ].set_xlim( 0.0, 1.0 )
        axDict[ refCn ].set_ylim( minCn, maxCn + 1.0 )

def drawCnvPlot( sample, options ):
    sampleName = sample.attrib[ 'sampleName' ]
    #print sampleName
    options.out = os.path.join( options.outdir, 'cnv_%s' %sampleName  )
    fig, pdf = libplot.initImage( 11.0, 3.25, options )

    title = "Copy Number Variation between %s and the reference %s" % (sampleName, sample.attrib['referenceName'] )
    
    cnvDict, minCn, maxCn = getSampleData( sample )
    axDict = setAxes( fig, cnvDict.keys(), options )
    for r in axDict:
        if r != 'bg':
            drawOneCnvPlot( r, axDict[ r ], cnvDict[ r ], options, minCn, maxCn )
    drawAxisLabels( axDict, cnvDict, options, title, maxCn )
    setAxisLimits( axDict, minCn, maxCn )
    libplot.writeImage( fig, pdf, options )

def readfile( file, filteredSamples ):
    xmltree = ET.parse( file )
    root = xmltree.getroot()
    allsamples = root.findall( 'statsForSample' )
    samples = []
    for s in allsamples:
        name = s.attrib[ 'sampleName' ]
        if name != "ROOT" and name != "" and name not in filteredSamples:
            samples.append( s )
    return samples

def checkOptions( args, options, parser ):
    if len( args ) < 1:
        parser.error( "No input file provided!\n" )
    if not os.path.exists( args[0] ):
        parser.error( "Input file %s does not exist.\n"  %args[0] )
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split('-')
    else:
        options.filteredSamples = []

def initOptions( parser ):
    parser.add_option( '--outdir', dest='outdir', default='.', help='Output directory' )
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def main():
    usage = ('Usage: %prog [options] inputCopyNumberStats.xml')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    #Read in input file:
    samples = readfile( args[0], options.filteredSamples )
    for sample in samples:
        drawCnvPlot( sample, options )

if __name__ == "__main__":
    main()
