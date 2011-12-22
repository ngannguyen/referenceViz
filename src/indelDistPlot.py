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

#from numpy import *
from numpy import array, log
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import LogLocator
from matplotlib.font_manager import FontProperties

def getSample(samples, name):
    for s in samples:
        if s.attrib[ 'sampleName' ] == name:
            return s
    return None

def getLeaves( allsamples, filteredSamples ):
    #Return only samples that are leaves in the tree
    samples = []
    for sample in allsamples:
        name = sample.attrib[ 'sampleName' ]
        if name != "ROOT" and name != "" and name not in filteredSamples:
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
    axleft = 0.1
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

def getFreq( dist, proportion, culm ):
    freqDict = {}
    totalbases = sum(dist)
    if totalbases == 0:
        return [], []
    for v in dist:
        if v in freqDict.keys():
            freqDict[ v ] += 1
        else:
            freqDict[ v ] = 1
    
    xdata = sorted( freqDict )
    ydata = []
    for x in xdata:
        y = freqDict[ x ]
        if proportion:
            y *= x
            #y = 100.0*y/totalbases #Percentage
        ydata.append( y )
    
    if culm: #Cumulative data
        for i in xrange( len(ydata) -1, -1, -1 ):
            ydata[i] = sum(ydata[:i+1])
    
    #HACK:
    if proportion and not culm:
        index = len(xdata)
        for i, x in enumerate(xdata):
            if x > 100:
                index = i
                break
        xdata = xdata[:index]
        ydata = ydata[:index]
    #if xlogscale == "true":
    #    xdata = log( array(xdata) )
    #if ylogscale == "true":
    #    ydata = log( array(ydata) )

    return xdata, ydata

def getLargeIndelProp(xdata, ydata): 
    large = 0
    total = ydata[ len(ydata) - 1 ]
    index = -1
    for i, x in enumerate(xdata):
        if xdata[i] >= 1000:
            index = i -1
            break
    if index >= 0 and total > 0:
        return (total - ydata[index])*100.0/total
    return 0

def drawData( axesList, samples, samplesPerPlot, options, proportion, culm ):
    largeIns = [] #List of proportion of total indel bases that indels >= 1000bp take up, each element is for each sample
    largeDels = []

    if len(axesList) %2 != 0:
        sys.stderr.write( 'Number of axes must be even. Got %d\n' %len(axesList) )
        sys.exit( 1 )

    colors = libplot.getColors1()
    if len(samples) < 1:
        return
    if samples[0].attrib["referenceName"] == "reference":
        colors.pop(0)
    elif samples[0].attrib["referenceName"] == 'hg19':
        colors.pop(1)

    #styles = []

    c = -1
    textsize = 'x-small'
    linesDict = {}
    labelsDict = {}
    xmax = float('-inf')
    ymax = float('-inf')
    xmin = float('inf')
    ymin = float('inf')
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
            #insXdata, insYdata = getFreq( insDist, options.xlogscale, options.ylogscale )
            insXdata, insYdata = getFreq( insDist, proportion, culm )
            delDist = [int(val) for val in sample.attrib[ 'deletionSizeDistribution' ].split()]
            #delXdata, delYdata = getFreq( delDist, options.xlogscale, options.ylogscale )
            delXdata, delYdata = getFreq( delDist, proportion, culm )

            #LARGE INDELS, FOR paper STATS, not related to the plot:
            if proportion and culm:
                largeIns.append( getLargeIndelProp(insXdata, insYdata) )
                largeDels.append( getLargeIndelProp(delXdata, delYdata) )

            c += 1
            il = insAxes.plot( insXdata, insYdata, color=colors[c] )
            dl = delAxes.plot( delXdata, delYdata, color=colors[c] )

            inslines.append( il )
            dellines.append( dl )
            
            insXmax = xmax
            delXmax = xmax
            if len(insXdata) >0:
                insXmax = max(insXdata)
            if len(delXdata) > 0:
                delXmax = max(delXdata)
            xmax = max( [xmax, insXmax, delXmax] )

            insYmax = ymax
            delYmax = ymax
            if len(insYdata) >0:
                insYmax = max(insYdata)
            if len(delYdata) > 0:
                delYmax = max(delYdata)
            ymax = max( [ymax, insYmax, delYmax] )

            insXmin = xmin
            delXmin = xmin
            if len(insXdata) >0:
                insXmin = min(insXdata)
            if len(delXdata) > 0:
                delXmin = min(delXdata)
            xmin = min( [xmin, insXmin, delXmin] )

            insYmin = ymin
            delYmin = ymin
            if len(insYdata) >0:
                insYmin = min(insYdata)
            if len(delYdata) > 0:
                delYmin = min(delYdata)
            ymin = min( [ymin, insYmin, delYmin] )

            #xmax = max([xmax, max(insXdata), max(delXdata)])
            #ymax = max([ymax, max(insYdata), max(delYdata)])
        
        linesDict[ i ] = inslines
        labelsDict[ i ] = sampleNames
        linesDict[ i + len(axesList)/2 ] = dellines
        labelsDict[ i + len(axesList)/2 ] = sampleNames
        #fontp = FontProperties()
        #fontp.set_size( 'x-small' )
        if i == 0:
            insAxes.set_title( 'Insertions' )    
            delAxes.set_title( 'Deletions' )
        
    for i in range( len(axesList) ):
        axes = axesList[ i ]
        if options.xlogscale == "true":
            axes.set_xscale('log')

        #if options.ylogscale == "true" and not proportion:
        if options.ylogscale == "true":
            axes.set_yscale('log')
        libplot.editSpine( axes )
        
        axes.set_xlabel('Length (bp)', size = textsize)
        if not proportion:
            axes.set_ylabel('Event number', size = textsize)
        else:
            axes.set_ylabel('Number of positions', size = textsize)
        axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
        #if options.xlogscale == "true":
        #    axes.set_xlabel('Log 2 of length (bp)', size = textsize)
        #else:
        #    axes.set_xlabel('Length (bp)', size = textsize)
        #if options.ylogscale == "true":
        #    axes.set_ylabel('Log 2 of count', size = textsize)
        #else:
        #    axes.set_ylabel('Count', size = textsize)

        #Legend
        legend = axes.legend( linesDict[ i ], [ libplot.properName(n) for n in labelsDict[ i ]], 'upper right', ncol=3 )
        for t in legend.get_texts():
            t.set_fontsize('x-small')
        legend._drawFrame = False

        if options.xlogscale == "true":
            scale = len(str(xmax)) -1
            xticks = [ 10**x for x in range(scale + 1) ]
            axes.set_xticks( xticks )
        #if options.ylogscale == "true" and not proportion:
        if options.ylogscale == "true":
            scale = len(str(ymax)) -1
            yticks = [ 10**y for y in range(scale + 1) ]
            axes.set_yticks( yticks )


        for label in axes.get_xticklabels():
            #label.set_rotation(75)
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

        axes.set_ylim( ymin, ymax )
        if proportion and not culm:
            axes.set_xlim( xmin, 100 )
        else:
            axes.set_xlim( xmin, xmax )

    #PRINT THE LARGE INDEL STATS:
    if proportion and culm:
        sys.stderr.write("largeIndelStats\n")
        sys.stderr.write("Large insertions: %f\n" %( sum(largeIns)/len(largeIns) ))
        sys.stderr.write("Large deletions: %f\n" %( sum(largeDels)/len(largeDels) ))
        largeIndels = [ (largeIns[i] + largeDels[i])/2.0 for i in range(len(largeIns)) ]
        sys.stderr.write("IndelsAverage: %f\n" %( sum(largeIndels)/len(largeIndels) ))

    return

def drawPlots( options, samples, outname, proportion, culm ):
    #sort samples:
    #samples = sorted( samples, key=lambda s: s.attrib[ 'sampleName' ] )

    if len(samples) < 1:
        return

    refname = samples[0].attrib[ 'referenceName' ]
    sys.stderr.write("%s\n" %refname)
    #options.out = os.path.join( options.outdir, 'indelDist_' + refname )
    options.out = os.path.join( options.outdir, 'indelDist_' + outname )
    if proportion:
        options.out = os.path.join( options.outdir, 'indelDist2_' + outname )
    if culm:
        options.out += '_culm'
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    
    samplesPerPlot = 10
    axesList = setAxes( fig, len(samples), samplesPerPlot )

    lines = drawData( axesList, samples, samplesPerPlot, options, proportion, culm )

    libplot.writeImage( fig, pdf, options )


def readfiles( files ):
    statsList = [] #each element represents each xml file
    names = []
    for f in files:
        xmltree = ET.parse( f )
        statsList.append( xmltree )
        names.append( os.path.basename(f).lstrip('pathStats_').rstrip('.xml') )

    return statsList, names

def initOptions( parser ):
    #parser.add_option('--title', dest='title', default='',
    #                  help='')
    #parser.add_option('--ycutoff', dest='ycutoff', default=0.5, type='float',
    #                   help='Only points with y values from ycutoff to 1 are displayed')
    parser.add_option('--outdir', dest='outdir', default='.', help='Output directory')
    #parser.add_option('--samplesOrder', dest="samplesOrder", default="reference,hg19,apd,cox,dbb,mann,mcf,qbl,ssto,NA12891,NA12892,NA12878,NA19239,NA19238,NA19240,panTro2", help="Samples order")
    parser.add_option('--samplesOrder', dest="samplesOrder", default="reference,hg19,apd,cox,dbb,mann,mcf,qbl,ssto,venter,NA12892,NA12878,NA19239,NA19238,NA19240,nigerian,yanhuang,panTro3", help="Samples order")
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
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
    options.samplesOrder = options.samplesOrder.split(',')
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split('-')
    else:
        options.filteredSamples = []
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

    statsList, names = readfiles( options.files )

    #Filter out non-leaf samples:
    statsListF = []
    for xmltree in statsList:
        root = xmltree.getroot()
        allsamples = root.findall( 'statsForSample' )

        samples = getLeaves( allsamples, options.filteredSamples )
        statsListF.append( samples )

    for i, samples in enumerate(statsListF):
        outname = names[i]
        sortedSamples = []
        if len(options.samplesOrder) == 0:
            sortedSamples = sorted(samples, key=lambda s:s.name)
        else:
            for name in options.samplesOrder:
                for sample in samples:
                    if sample.attrib['sampleName'] == name:
                        sortedSamples.append( sample )
        #print outname
        drawPlots( options, sortedSamples, outname, False, False )
        drawPlots( options, sortedSamples, outname, True, False )
        drawPlots( options, sortedSamples, outname, True, True )

    #if len( statsListF ) >= 2:
    #    for i in range( len(statsListF) -1 ):
    #        for j in range( i + 1, len(statsListF) ):
    #            drawCompareN50Plot( options, statsListF[i], statsListF[j] )


if __name__ == '__main__':
    main()

