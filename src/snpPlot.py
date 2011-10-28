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
import matplotlib.pylab as pylab
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

class Sample():
    def __init__( self, xmlnode ):
        self.name = xmlnode.attrib[ 'sampleName' ]
        self.refname = xmlnode.attrib[ 'referenceName' ]
        self.totalErrors = int( xmlnode.attrib[ 'totalErrors' ] )
        self.totalCalls = int( xmlnode.attrib[ 'totalCalls' ] )
        self.errPerSite = 0
        if self.totalCalls != 0:
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
            if name != 'aggregate' and name != 'ROOT' and name != '' and name not in options.filteredSamples:
                samples.append( Sample( s ) )
        statsList.append( samples )
    return statsList

def initOptions( parser ):
    parser.add_option( '--outdir', dest='outdir', default='.', help='Output directory' ) 
    parser.add_option( '--numOutliners', dest='numOutliners', default=1, help='Number of outliners' ) 
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def checkOptions( args, options, parser ):
    if len( args ) < 2:
        parser.error( 'Please provide two snpStats xml files.\n' )
    options.files = []
    for f in args:
        if not os.path.exists(f):
            parser.error( 'File %s does not exist.\n'  % f )
        options.files.append( f )
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split('-')
    else:
        options.filteredSamples = []

def setAxes(fig, range1, range2):
    axleft = 0.12
    axright = 0.95
    axwidth = axright - axleft
    axbottom = 0.15
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.015

    h1 = (axheight - margin)*(range1/(range1 + range2))
    h2 = axheight - margin - h1

    ax2 = fig.add_axes( [axleft, axbottom, axwidth, h2] )
    ax = fig.add_axes( [axleft, axbottom + h2 + margin, axwidth, h1] )
    return ax, ax2

def drawSnpPlot( options, samples1, samples2 ):
    #Sorted in decreasing order of errorPerSite in samples1
    samples1 = sorted( samples1, key=lambda s:s.errPerSite, reverse=True )
    if len( samples1 ) < 1:
        return
    
    refname1 = samples1[0].refname
    refname2 = samples2[0].refname

    y1data = [ s.errPerSite for s in samples1 ]
    
    xticklabels = [ s.name for s in samples1 ]
    y2data = []
    for name in xticklabels:
        for s in samples2:
            if s.name == name or (name == refname2 and s.name == refname1):
                y2data.append(s.errPerSite)
                break
    if len(xticklabels) != len(y2data):
        sys.stderr.write("Input file 1 and 2 do not have the same set of samples\n")
        sys.exit( 1 )

    for i in range(len(xticklabels)):
        if xticklabels[i] == refname2:
            xticklabels[i] = "ref/hg19"
    
    #add the average column:
    num = options.numOutliners
    if len(y1data) > num:
        y1avr = sum(y1data[num:])/float(len(y1data) - num)
        y1data.append(y1avr)
        xticklabels.append('average')
        y2avr = sum(y2data[num:])/float(len(y2data) - num)
        y2data.append(y2avr)

    minOutlier = min( [min(y2data[0:num]), min(y1data[0:num])] ) - 0.001 
    maxOutlier = max( [max(y2data[0:num]), max(y1data[0:num])] ) + 0.001
    minMajority = min( [min(y2data[num:]), min(y1data[num:])] ) - 0.001
    maxMajority = max( [max(y2data[num:]), max(y1data[num:])] ) + 0.001

    basename = os.path.basename(options.files[0])
    #options.out = os.path.join( options.outdir, '%s' %(basename.split('_')[0]) )
    options.out = os.path.join( options.outdir, '%s' %( basename.lstrip('snpStats').lstrip('_').rstrip('.xml') ) )
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    ax, ax2 = setAxes(fig, maxOutlier - minOutlier, maxMajority - minMajority)

    #Plot the outliers:
    l2 = ax.plot( y2data, marker='.', markersize=14.0, linestyle='none', color="#1F78B4" )
    l1 = ax.plot( y1data, marker='.', markersize=14.0, linestyle='none', color="#E31A1C" )

    ax2.plot( y2data, marker='.', markersize=14.0, linestyle='none', color="#1F78B4" )
    ax2.plot( y1data, marker='.', markersize=14.0, linestyle='none', color="#E31A1C" )
  
    #Legend
    fontP = FontProperties()
    fontP.set_size("x-small")
    legend = ax.legend([l1, l2], [libplot.properName(refname1), libplot.properName(refname2)], 'upper right', numpoints=1, prop=fontP)
    legend._drawFrame = False

    d = .0001 # how big to make the diagonal lines in axes coordinates
    ax.plot( (-1, 0), (minOutlier +d, minOutlier - d), color = "k", clip_on=False )
    ax2.plot( (-1, 0), (maxMajority +d, maxMajority - d), color = "k", clip_on=False )
    
    ax.set_ylim( minOutlier, maxOutlier ) # outliers only
    ax.set_xlim( -0.5, len(xticklabels) -0.5 )
    dummyxticklabels = [ "" for l in xticklabels ]
    ax.set_xticklabels( dummyxticklabels )
    
    #ytickpositions = ax.yaxis.get_majorticklocs()
    #ytickpositions = [ ytickpositions[i] for i in range(0,len(ytickpositions),2) ]
    #Make sure the y ticks of the top plot (the outlier plot) is the same with the other plot:
    step = 0.001
    ytickpositions = []
    #ytickpos = 0.013
    ytickpos = 0
    while ytickpos < maxOutlier:
        if ytickpos >= minOutlier:
            ytickpositions.append(ytickpos)
        ytickpos += step
    ax.set_yticks( ytickpositions )
        
    #ax.xaxis.set_major_locator( NullLocator() )
    #ax.xaxis.set_major_formatter( NullFormatter() )
    
    ax2.set_ylim( minMajority, maxMajority )
    ax2.set_xlim( -0.5, len(xticklabels) -0.5 )

    # hide the spines between ax and ax2
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position( 'left' )
    ax.xaxis.set_ticks_position( 'none' )

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.yaxis.set_ticks_position( 'left' )

    ax2.set_xticks( range( 0, len(xticklabels) ) )
    properxticklabels = [ libplot.properName(l) for l in xticklabels ]
    ax2.set_xticklabels( properxticklabels )
    #Make sure the x ticks of the top plot is the same with the other plot:
    ax.set_xticks( range(0, len(xticklabels)) )

    for label in ax2.xaxis.get_ticklabels():
        label.set_rotation( 90 )
   
    ax.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax2.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax2.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)

    ax2.set_xlabel( 'Samples' )
    ax2.set_ylabel( 'Snps per site' )
    title = 'SNPs'
    ax.set_title( title )
    
    libplot.writeImage( fig, pdf, options )
   
def getSample( samples, name ):
    for s in samples:
        if s.name == name:
            return s
    return None

def main():
    usage = ('Usage: %prog [options] file1.xml file2.xml\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    statsList = readfiles( options )

    drawSnpPlot( options, statsList[0], statsList[1] )
    

if __name__ == "__main__":
    main()
