#!/usr/bin/env python

"""
Non-linear breakpoints plot (xaxis = samples, yaxis = number of non-linear bps in that sample
nknguyen at soe dot ucsc dot edu
Input: pathStats.xml files
"""
import os, sys
from optparse import OptionParser
import xml.etree.ElementTree as ET

import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.font_manager import FontProperties

def drawPlot(samplesList, sampleNames, options):
    options.out = os.path.join(options.outdir, "nonLinearBp")
    fig, pdf = libplot.initImage(12.0, 8.0, options)
    axes = fig.add_axes([0.09, 0.2, 0.9, 0.6])
   
    list1 = samplesList[0]
    list2 = samplesList[1]
    if len(list1) < 1 or len(list2) < 1:
        return
    refname1 = list1[0].attrib['referenceName']
    refname2 = list2[0].attrib['referenceName']

    lines = []

    barwidth = 0.3
    y1data = []
    y2data = []
    for sample in sampleNames:
        for s in list1:
            if sample == s.attrib['sampleName']:
                y1data.append( int(s.attrib['totalIntraJoin']) )
        for s in list2:
            if sample == s.attrib['sampleName']:
                y2data.append( int(s.attrib['totalIntraJoin']) )
    x1data = range( len(y1data) )
    x2data = [ x+ barwidth for x in x1data]
     
    colors =["#1F78B4", "#E31A1C"]
    l1 = axes.bar( x1data, y1data, barwidth, color = colors[0], ec='w')
    lines.append(l1[0])
    l2 = axes.bar( x2data, y2data, barwidth, color = colors[1], ec='w')
    lines.append(l2[0])

    libplot.editSpine(axes)
    axes.set_title("Non-linear Breakpoints")
    
    #set ticks
    xlabels = [ libplot.properName(name) for name in sampleNames ]
    fontP = FontProperties()
    fontP.set_size('small')
    pyplot.xticks(x2data, xlabels, rotation=45, fontproperties=fontP)
    pyplot.yticks( fontproperties = fontP )
    pyplot.xlabel("Samples")
    pyplot.ylabel("Number of breakpoints")
    axes.xaxis.set_ticks_position('bottom')
    axes.yaxis.set_ticks_position('left')
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)
    legend = axes.legend( lines, [libplot.properName(refname1), libplot.properName(refname2)], prop=fontP, loc="best" )
    legend._drawFrame = False

    libplot.writeImage(fig, pdf, options)

def readfiles(files, filteredSamples):
    samplesList = []
    for f in files:
        tree = ET.parse(f)
        root = tree.getroot()
        samples = []
        for sample in root.findall( 'statsForSample' ):
            name = sample.attrib['sampleName']
            if name != "ROOT" and name != "" and name not in filteredSamples:
                samples.append(sample)
        samplesList.append( samples )
    return samplesList

def initOptions( parser ):
    parser.add_option('--title', dest='title', default='Coverage statistics', help='Based title of the plots, default=%default')
    parser.add_option('--ycutoff', dest='ycutoff', default=0.9, type='float')
    parser.add_option('--outdir', dest='outdir', default='.')
    parser.add_option('--samples', dest='samples', default='apd,cox,dbb,mann,mcf,qbl,ssto,venter,nigerian,yanhuang,NA12892,NA12878,NA19239,NA19238,NA19240,panTro3', help="Comma separated list of samples in the desired order")
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def checkOptions( args, options, parser ):
    if len(args) < 2:
        parser.error('Please specify two pathStats*.xml files\n')
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
    usage = ('Usage: %prog [options] pathStats1.xml pathStats2.xml\n\n'
             '%prog takes in pathStats*.xml files and create the non-linear breakpoint dist. plots\n')
    parser = OptionParser(usage = usage)
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()

    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )
    
    samplesList = readfiles( options.files, options.filteredSamples )
    items = options.samples.split(',')
    sampleNames = []
    for item in items:
        if item not in options.filteredSamples:
            sampleNames.append(item)
    list1 = samplesList[0]
    list2 = samplesList[1]
    drawPlot(samplesList, sampleNames, options)



if __name__ == "__main__":
    main()
