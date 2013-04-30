#!/usr/bin/env python2.6

"""
Create Snp plots
nknguyen at soe dot ucsc dot edu
Jun 15 2011
"""
import os, sys, re
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
        self.dels = float(xmlnode.attrib['totalDeletionPerAlignedBase'])
        self.ins = float(xmlnode.attrib['totalInsertionPerAlignedBase'])
        self.indel = float(xmlnode.attrib['totalInsertionAndDeletionPerAlignedBase'])
        #add cases where it's both insertion & deletion (indel) to insertion & deletion:
        self.dels += self.indel
        self.ins += self.indel

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
    #parser.add_option( '--numOutliners', dest='numOutliners', type='int', default=0, help='Number of outliners' ) 
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

def drawPlot_noOutgroup( options, samples1, samples2, type ):
    #Sorted in decreasing order of errorPerSite in samples1
    if type == 'insertion':
        samples1 = sorted( samples1, key=lambda s:s.ins, reverse=True )
    else:
        samples1 = sorted( samples1, key=lambda s:s.dels, reverse=True )
    if len( samples1 ) < 1:
        return
    
    refname1 = samples1[0].refname
    refname2 = samples2[0].refname

    y1data = [ s.ins for s in samples1 ]
    if type == 'deletion':
        y1data = [ s.dels for s in samples1 ]
    xticklabels = [ s.name for s in samples1 ]
    
    #indel of refname1 w.r.t itself (0)
    y1data.append(0)
    xticklabels.append(refname1)

    y2data = []
    for name in xticklabels:
        if name == refname2:#indel of refname2 w.r.t itself (0)
            y2data.append(0)
        for s in samples2:
            if s.name == name:
                if type == 'insertion':
                    y2data.append(s.ins)
                else:
                    y2data.append(s.dels)
                break
    
    if len(xticklabels) != len(y2data):
        sys.stderr.write("Input file 1 and 2 do not have the same set of samples\n")
        sys.exit( 1 )

    #add the average column:
    num = 1
    y1avr = sum(y1data)/float(len(y1data) - 1)
    y1data.append(y1avr)
    xticklabels.append('average')
    y2avr = sum(y2data)/float(len(y2data) - 1)
    y2data.append(y2avr)
    print "%s Average: %s %f, %s %f" %(type, refname1, y1avr, refname2, y2avr)

    minMajority = min( [min(y2data), min(y1data)] ) - 0.0001
    maxMajority = max( [max(y2data), max(y1data)] ) + 0.0001

    basename = os.path.basename(options.files[0])
    options.out = os.path.join( options.outdir, '%s_%s' %( type, basename.lstrip('pathStats').lstrip('_').rstrip('.xml') ) )
    fig, pdf = libplot.initImage( 11.2, 10.0, options )
    ax2 = fig.add_axes( [0.15, 0.15, 0.8, 0.8] )

    l2 = ax2.plot( y2data, marker='.', markersize=14.0, linestyle='none', color="#E31A1C" )
    l1 = ax2.plot( y1data, marker='.', markersize=14.0, linestyle='none', color="#1F78B4" )
    
    #Legend
    fontP = FontProperties()
    fontP.set_size("x-small")
    legend = ax2.legend([l1, l2], [libplot.properName(refname1), libplot.properName(refname2)], 'upper right', numpoints=1, prop=fontP)
    legend._drawFrame = False
            
    ax2.set_ylim( minMajority, maxMajority )
    ax2.set_xlim( -0.5, len(xticklabels) -0.5 )

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.yaxis.set_ticks_position( 'left' )

    ax2.set_xticks( range( 0, len(xticklabels) ) )
    properxticklabels = [ libplot.properName(l) for l in xticklabels ]
    ax2.set_xticklabels( properxticklabels )

    for label in ax2.xaxis.get_ticklabels():
        label.set_rotation( 90 )
   
    ax2.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax2.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)

    ax2.set_xlabel( 'Samples' )
    title = 'Deletions'
    #if type == 'insertion':
    if type == 'insertion':
        ax2.set_ylabel( 'Insertions per site' )
        title = 'Insertions'
    else:
        ax2.set_ylabel( 'Deletions per site' )
    ax2.set_title( title )
    
    libplot.writeImage( fig, pdf, options )
   
def drawPlot( options, samples1, samples2, type ):
    #Sorted in decreasing order of errorPerSite in samples1
    if type == 'insertion':
        samples1 = sorted( samples1, key=lambda s:s.ins, reverse=True )
    else:
        samples1 = sorted( samples1, key=lambda s:s.dels, reverse=True )
    if len( samples1 ) < 1:
        return
    
    #remove outgroupSample:
    outgroupSample = options.outgroup
    for i, s in enumerate(samples1):
        if re.search(options.outgroup, s.name):
            outgroupSample = samples1.pop(i)
            break

    refname1 = samples1[0].refname
    refname2 = samples2[0].refname

    y1data = [ s.ins for s in samples1 ]
    if type == 'deletion':
        y1data = [ s.dels for s in samples1 ]
    xticklabels = [ s.name for s in samples1 ]
    
    #indel of refname1 w.r.t itself (0)
    y1data.append(0)
    xticklabels.append(refname1)

    y2data = []
    for name in xticklabels:
        if name == refname2:#indel of refname2 w.r.t itself (0)
            y2data.append(0)
        for s in samples2:
            if s.name == name:
                if type == 'insertion':
                    y2data.append(s.ins)
                else:
                    y2data.append(s.dels)
                break
    
    if len(xticklabels) != len(y2data):
        sys.stderr.write("Input file 1 and 2 do not have the same set of samples\n")
        sys.exit( 1 )

    #add the average column:
    num = 1
    y1avr = sum(y1data)/float(len(y1data) - 1)
    y1data.append(y1avr)
    xticklabels.append('average')
    y2avr = sum(y2data)/float(len(y2data) - 1)
    y2data.append(y2avr)
    print "%s Average: %s %f, %s %f" %(type, refname1, y1avr, refname2, y2avr)

    #Add outgroup:
    if outgroupSample:
        samples1.append(outgroupSample)
        if type == 'insertion':
            y1data.append( outgroupSample.ins )
        else:
            y1data.append( outgroupSample.dels )
        for s in samples2:
            if re.search(options.outgroup, s.name):
                if type == 'insertion':
                    y2data.append(s.ins)
                else:
                    y2data.append(s.dels)
        xticklabels.append(options.outgroup)

    minMajority = min( [min(y2data), min(y1data)] ) - 0.0001
    maxMajority = max( [max(y2data), max(y1data)] ) + 0.0001

    basename = os.path.basename(options.files[0])
    options.out = os.path.join( options.outdir, '%s_%s' %( type, basename.lstrip('pathStats').lstrip('_').rstrip('.xml') ) )
    fig, pdf = libplot.initImage( 11.2, 10.0, options )
    #ax, ax2 = setAxes(fig, maxOutlier - minOutlier, maxMajority - minMajority)
    ax2 = fig.add_axes( [0.15, 0.15, 0.8, 0.8] )

    l2 = ax2.plot( y2data, marker='.', markersize=14.0, linestyle='none', color="#E31A1C" )
    l1 = ax2.plot( y1data, marker='.', markersize=14.0, linestyle='none', color="#1F78B4" )
    
    #Legend
    fontP = FontProperties()
    fontP.set_size("x-small")
    legend = ax2.legend([l1, l2], [libplot.properName(refname1), libplot.properName(refname2)], 'upper right', numpoints=1, prop=fontP)
    legend._drawFrame = False
            
    ax2.set_ylim( minMajority, maxMajority )
    ax2.set_xlim( -0.5, len(xticklabels) -0.5 )

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.yaxis.set_ticks_position( 'left' )

    ax2.set_xticks( range( 0, len(xticklabels) ) )
    properxticklabels = [ libplot.properName(l) for l in xticklabels ]
    ax2.set_xticklabels( properxticklabels )

    for label in ax2.xaxis.get_ticklabels():
        label.set_rotation( 90 )
   
    ax2.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax2.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)

    ax2.set_xlabel( 'Samples' )
    title = 'Deletions'
    #if type == 'insertion':
    if type == 'insertion':
        ax2.set_ylabel( 'Insertions per site' )
        title = 'Insertions'
    else:
        ax2.set_ylabel( 'Deletions per site' )
    ax2.set_title( title )
    
    libplot.writeImage( fig, pdf, options )
   
def main():
    usage = ('Usage: %prog [options] file1.xml file2.xml\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )
    parser.add_option('--outgroup', dest='outgroup', help='Name of outgroup sample')

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    statsList = readfiles( options )

    #drawPlot( options, statsList[0], statsList[1], 'insertion' )
    #drawPlot( options, statsList[0], statsList[1], 'deletion' )
    drawPlot_noOutgroup( options, statsList[0], statsList[1], 'insertion' )
    drawPlot_noOutgroup( options, statsList[0], statsList[1], 'deletion' )
    

if __name__ == "__main__":
    main()
