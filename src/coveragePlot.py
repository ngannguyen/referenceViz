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
    axes.set_title("Sample Coverage") #TO BE NAMED!!!
    pyplot.xlabel("Samples")
    if isAbs:
        pyplot.ylabel("Number of positions")
    else:
        pyplot.ylabel("Proportion of total positions")

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

#================= Latex table of mapped bases of each sample onto the reference vs hg19 ===========
def tabHeader(f, ref1, ref2):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|r|r|r|r}\n")
    #f.write("\\multicolumn{6}{c}{%s} \\\\\n" %title)
    #f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("Sample & Repeat & %s & %s & Total\\\\\n" %(libplot.properName(ref1), libplot.properName(ref2)))
    #f.write("Sample & Reads & Total Bases & \\%%Repeats & SNP Rate & Overal Snp Rate\\\\\n")
    f.write("\\hline\n")

def tab(f, stats, sample2repeat):
    altColor = 1 
    for sample in stats:
        repeat = 'NA'
        repeatPc = ''
        if sample.name in sample2repeat:
            repeat = libplot.prettyInt( sample2repeat[sample.name][1] )
            repeatPc = "(%.2f \\%%)" %sample2repeat[sample.name][2]
        otherRef = libplot.prettyInt(sample.otherReferenceBasesMapped)
        otherRefPc = "%.2f" % (100.0*sample.otherReferenceBasesMapped/sample.totalBases)
        ref = libplot.prettyInt(sample.referenceBasesMapped)
        refPc = "%.2f" % (100.0*sample.referenceBasesMapped/sample.totalBases)
        total = libplot.prettyInt(sample.totalBases)
        sampleName = libplot.properName(sample.name)
        if altColor == 1:
            f.write("%s & %s %s & %s (%s \\%%) & %s (%s \\%%) & %s \\\\\n" %(sampleName, repeat, repeatPc, otherRef, otherRefPc, ref, refPc, total ))
        else:
            f.write("\\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s %s & \\cellcolor[gray]{0.9} %s (%s \\%%) & \\cellcolor[gray]{0.9} %s (%s \\%%) & \\cellcolor[gray]{0.9} %s \\\\\n" %(sampleName, repeat, repeatPc, otherRef, otherRefPc, ref, refPc, total ))
        altColor = 1 - altColor
    f.write("\\hline\n")

def drawCompareCoverageTab( options, stats, sample2repeat ):
    prefix = "cmpCoverageTab"
    outfile = os.path.join( options.outdir, "%s%s_%s.tex" %(prefix, stats[0].referenceName, stats[0].otherReferenceName) )
    f = open(outfile, 'w')
    libplot.writeDocumentStart(f)
    tabHeader(f, stats[0].otherReferenceName, stats[0].referenceName)
    tab(f, stats, sample2repeat)
    captionStr = ""
    label = ""
    libplot.tableCloser(f, captionStr, label)
    libplot.writeDocumentEnd(f)
    f.close()

#================= Mapped bases of each sample onto the reference vs hg19 ===========
def drawCompareData( axes, options, stats, isAbs ):
    if len(stats) == 0:
        return
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = [ stats[0].otherReferenceName, stats[0].referenceName, "total" ]

    barwidth = 0.25
    #X data:
    x3data = []
    #avgIndex = -1

    currx = -1
    #xVer = [] #location (x) of vertical lines to separate between human samples | avr, all | chimp
    for i,s in enumerate( stats ):
        #if s.name == 'average':
        #    avgIndex = i
        if s.name == 'average' or s.name == 'panTro3':
            currx += 1 + 1.5*barwidth
            #xVer.append( currx - (1.0 + 1.5*barwidth - 3*barwidth)/2.0 )
        else:
            currx += 1
        x3data.append( currx )
    
    #print x1data
    x2data = [ x + barwidth for x in x3data ]
    x1data = [ x + barwidth for x in x2data ]

    if isAbs:
        y1data = [ sample.otherReferenceBasesMapped for sample in stats ]
        y2data = [ sample.referenceBasesMapped for sample in stats ]
        y3data = [ sample.totalBases for sample in stats ]
    else:
        y1data = [ 100.0*sample.otherReferenceBasesMapped/sample.totalBases for sample in stats ]
        y2data = [ 100.0*sample.referenceBasesMapped/sample.totalBases for sample in stats ]
        y3data = [ 100.0*sample.totalBases/sample.totalBases for sample in stats ]
    
    #Average aggregate data:
    #if avgIndex > 0:
    #    y1data[ avgIndex ] /= float(avgIndex)
    #    y2data[ avgIndex ] /= float(avgIndex)
    #    y3data[ avgIndex ] /= float(avgIndex)

    colors =["#1F78B4", "#E31A1C", "#4DAF4A"] 
    #colors =["#1B9E77", "#D95F02", "#7570B3"] 
    #colors =["#EDF8B1", "#7FCDBB", "#2C7FB8"]
    #colors =["#A1DAB4", "#41B6C4", "#225EA8"]
    l1 = axes.bar( x1data, y1data, barwidth, color = colors[0], ec="w" ) 
    lines.append( l1[0] )
    l2 = axes.bar( x2data, y2data, barwidth, color = colors[1], ec="w" ) 
    lines.append( l2[0] )
    l3 = axes.bar( x3data, y3data, barwidth, color = colors[2], ec="w" )
    lines.append( l3[0] )

    libplot.editSpine( axes )
    axes.set_title("Sample Coverage") #TO BE NAMED

    #set ticks:
    samples = []
    for sample in stats:
        samples.append( libplot.properName(sample.name) )
    fontP = FontProperties()
    fontP.set_size('small')
    #pyplot.xticks( x + barwidth/2., samples, rotation=45, fontproperties=fontP )
    pyplot.xticks( x2data, samples, rotation=45, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )

    #HACK:
    yticks = range(2000000, 6000000, 500000)
    yticklabels = [ float(y)/1000000 for y in yticks ]
    axes.set_yticks(yticks)
    axes.set_yticklabels(yticklabels)
    
    pyplot.xlabel("Samples")
    pyplot.ylabel("Number of positions (in millions)")
    
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    miny = min( [min(y1data), min(y2data), min(y3data)] )
    miny = miny*0.9
    maxy = max([max(y1data), max(y2data), max(y3data)])
    
    #Draw vertical lines:
    #for x in xVer:
    #    axes.plot([x, x], [miny, maxy], color="#A8A8A8")

    axes.set_ylim( miny, maxy )
    axes.set_xlim(-0.5, max(x1data) + 0.5 )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)
    
    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.95, box.height*0.9] )

    legend = axes.legend( lines, [libplot.properName(n) for n in linenames], prop=fontP, loc="best", bbox_to_anchor=(0.2, 1) )
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

#================== DRAW ONLY BASES OF EACH SAMPLE THAT MAPPED TO CACTUS REF ==================
def drawCompareData2( axes, options, stats, isAbs ):
    if len(stats) == 0:
        return
    #if isAbs, draw absolute values. If not, draw proportion (relative values)
    lines = []
    linenames = [ stats[0].otherReferenceName, stats[0].referenceName, "total" ]

    #X data:
    x1data = []

    currx = -1
    for i,s in enumerate( stats ):
        if s.name == 'all':
            continue
        if s.name == 'average' or s.name == 'panTro3':
            currx += 1.5
        else:
            currx += 1
        x1data.append( currx )

    y1data = []
    for sample in stats:
        if sample.name == 'all':
            continue
        if isAbs:
            y1data.append( sample.referenceBasesMapped )
        else:
            y1data.append( 100.0*sample.referenceBasesMapped/sample.totalBases )
    
    barwidth = 0.6
    #barwidth = 0.25
    l1 = axes.bar( x1data, y1data, barwidth, color = "#E31A1C", ec="w" ) 
    lines.append( l1[0] )

    libplot.editSpine( axes )
    axes.set_title("Sample Coverage") #TO BE NAMED

    #set ticks:
    samples = []
    for sample in stats:
        if sample.name == 'all':
            continue
        samples.append( libplot.properName(sample.name) )
    fontP = FontProperties()
    fontP.set_size('small')
    pyplot.xticks( [x + barwidth/2.0 for x in x1data], samples, rotation=45, fontproperties=fontP )
    pyplot.yticks( fontproperties=fontP )

    #HACK:
    yticks = range(2000000, 6000000, 500000)
    yticklabels = [ float(y)/1000000 for y in yticks ]
    axes.set_yticks(yticks)
    axes.set_yticklabels(yticklabels)
    
    pyplot.xlabel("Samples")
    pyplot.ylabel("Number of positions (in millions)")
    
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )

    miny = min( y1data )
    miny = miny*0.9
    axes.set_ylim( miny, max(y1data) )
    axes.set_xlim(-0.5, max(x1data) + 0.5 )
    
    axes.yaxis.grid(b=True, color="#A8A8A8", linestyle='-', linewidth=0.25)

    #Legend:
    box = axes.get_position()
    axes.set_position( [box.x0, box.y0, box.width*0.95, box.height*0.9] )

    #legend = axes.legend( lines, [libplot.properName(n) for n in linenames], prop=fontP, loc="best", bbox_to_anchor=(0.2, 1) )
    #legend._drawFrame=False
    
    return 

def drawCompareCoveragePlot2( options, stats, isAbs ):
    if len(stats) == 0:
        return
    prefix = "cmpCoverage2_"
    if not isAbs:
        prefix = "cmpRelCoverage2_"
    options.out = os.path.join( options.outdir, "%s%s_%s" %(prefix, stats[0].referenceName, stats[0].otherReferenceName) )
    fig, pdf = libplot.initImage( 12.0, 8.0, options )
    axes = fig.add_axes( [0.09, 0.2, 0.9, 0.6] )

    drawCompareData2( axes, options, stats, isAbs )
    libplot.writeImage( fig, pdf, options )

#================ display data in scatter plot form, xaxis = number of sample covered, yaxis = number of bases
def drawScatter( axes, options, stats, type, cumulative ):
    if len(stats) < 4:
        return
    
    title = "Distribution of Positions Shared Among Samples"
    if cumulative:
        title = "Cumulative Distribution of Positions Shared Among Samples"
    axes.set_title(title) #TO BE NAMED
    
    #samples = ["panTro3", "minusOtherReference", "average", "reference", "hg19"]
    samples = ["reference", "hg19", "panTro3", "average"]
    
    if type == 'noHg19':
        samples = ["minusOtherReference"]
    xdata = range( 0, len(stats) -4 )
    #print xdata
    ydataList = []
    miny = float('inf')
    maxy = float('-inf')
    for name in samples:
        for s in stats:
            if s.name == name:
                ydata = s.baseCoverages[: len(stats) -4]
                if cumulative:
                    ydata = [ sum(ydata[i:]) for i in xrange( len(ydata) ) ]

                ydataList.append( ydata )
                miny = min( [miny, min(ydata)] )
                maxy = max( [maxy, max(ydata)] )
                break

    lines = []
    #colors = libplot.getColors0()
    colors =["#E31A1C", "#1F78B4", "#3D3D3D", "#4DAF4A"] #ConsensusRef, GRCh37, chimp, average
    c = -1
    offset = 0.12
    axes.set_yscale('log')
    #if type == 'noHg19':
    #    axes.set_yscale('log')

    for i in xrange( len(samples) ):
        xdatai = [x + offset*i for x in xdata]
        ydata = ydataList[i]
        c += 1
        if i == 0:
            axes.plot(xdatai[1:], ydata[1:], color="#CCCCCC", linestyle='-', linewidth=0.002)
        else:
            axes.plot(xdatai, ydata, color="#CCCCCC", linestyle='-', linewidth=0.002)
        l = axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=12.0, linestyle='none')
        lines.append(l)
    
    fontP = FontProperties()
    fontP.set_size('x-small')

    yrange = maxy - miny
    miny = miny - 10
    maxy = maxy + yrange*0.1
    
    xmin = -0.4
    xmax = len(stats) - 4 -1 + offset*len(samples) + offset
    libplot.editSpine(axes)
    
    axes.set_xticks( [ i + offset*(len(samples)/2.0 ) for i in range(0, len(stats) -4)] )
    axes.set_xticklabels( range(1, len(stats) -2) )
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    scale = len(str( int(maxy) )) - 1
    ylabel = "Number of positions"
    if type == "noHg19": 
        yticks = [ 10**y for y in range(scale + 1) ]
    else:
        #yticks = [ 10**y for y in range(scale + 2) ]
        yticks = []
        for y in range(scale + 1):
            for y1 in range(1,10):
                yticks.append(y1*(10**y))
    axes.set_yticks( yticks )
    minorLocator = LogLocator( base=10, subs = range(1, 10) )
    axes.yaxis.set_minor_locator( minorLocator )
    #else:
    #    yticks = range(0, int(maxy), 10**scale)
    #    yticklabels = [ y/(10**scale) for y in yticks ]
    #    axes.set_yticks( yticks )
    #    axes.set_yticklabels( yticklabels )
    #    ylabel += " (x%s)" %( libplot.prettyInt(10**scale) )
        #ylabel += " (in millions)"
    axes.set_xlim(xmin, xmax)
    if type == "noHg19":
        axes.set_ylim(miny, maxy)
    else:
        axes.set_ylim(10000, 1000000)#HACK

    if type != 'noHg19':
        legend = pyplot.legend( lines, [libplot.properName(s) for s in samples], numpoints=1, loc='lower right', prop=fontP )
        legend._drawFrame = False

    axes.set_xlabel( 'Number of samples' )
    #if type == "noHg19":
    #    ylabel += " (x %d)" %(10**(scale -1))
    axes.set_ylabel( ylabel )
    #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    return

def drawScatterPlot( options, stats, type, cumulative ):
    prefix = "coverageScatter_%s" %type
    if cumulative:
        prefix += "_culm"
    options.out = os.path.join( options.outdir, "%s" %(prefix) )
    fig, pdf = libplot.initImage( 12.0, 8.0, options )
    axes = fig.add_axes( [0.1, 0.15, 0.85, 0.75] )
    drawScatter( axes, options, stats, type, cumulative )
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

def readRepeatInfo(file):
    f = open(file, 'r')
    sample2repeat = {} #key = sample, vals = [totalBases, repetitiveBases, Percentage]
    f.readline()
    for line in f.readlines():
        items = line.strip().split('\t')
        if len(items) < 4:
            continue
        sample2repeat[ items[0] ] = [ int(items[1]), int(items[2]), float(items[3]) ]
    #Calculate average:
    avr = [0, 0, 0]
    numSamples = 0
    for sample in sample2repeat:
        if sample == 'panTro3':
            continue
        numSamples += 1
        for i, val in enumerate( sample2repeat[sample] ):
            avr[i] += val
    if numSamples > 0:
        for i in xrange( len(avr) ):
            avr[i] /= numSamples
    sample2repeat['average'] = avr
    return sample2repeat

def initOptions( parser ):
    parser.add_option('--title', dest='title', default='Coverage statistics', help='Based title of the plots, default=%default')
    parser.add_option('--ycutoff', dest='ycutoff', default=0.9, type='float')
    parser.add_option('--outdir', dest='outdir', default='.')
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')
    parser.add_option('--repeatInfo', dest='repeatInfo', default='/hive/users/nknguyen/reconGit/referenceScripts/data/seqRepeatPercentage.txt')
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
    sample2repeat = readRepeatInfo( options.repeatInfo )
    
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
        sortedstats = sorted( sortedstats, key=lambda s: s.totalBases, reverse=True )
        for k in specialcases:
            s = specialcases[k]
            if s:
                sortedstats.append( s )
        if len(sortedstats) > 0:
            drawCompareCoveragePlot( options, sortedstats, True )
            drawCompareCoverageTab( options, sortedstats, sample2repeat )

        #sort by totalBases, but include 'otherReference' sample
        sortedstats = []
        for i in xrange( len(stats) ):
            if stats[i].name in [ stats[i].referenceName, 'minusOtherReference' ]:
                continue
            if stats[i].name in specialcases:
                specialcases[ stats[i].name ] = stats[i]
            else:
                sortedstats.append( stats[i] )
        sortedstats = sorted( sortedstats, key=lambda s: s.totalBases, reverse=True )
        for k in specialcases:
            s = specialcases[k]
            if s:
                sortedstats.append( s )
        if len(sortedstats) > 0:
            drawCompareCoveragePlot2( options, sortedstats, True )
        
        #scatter plot
        drawScatterPlot( options, stats, 'consensusVShg19', False )
        drawScatterPlot( options, stats, 'consensusVShg19', True )
        drawScatterPlot( options, stats, 'noHg19', False )
        drawScatterPlot( options, stats, 'noHg19', True )


if __name__ == "__main__":
    main()





