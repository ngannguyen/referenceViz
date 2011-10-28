#!/usr/bin/env python

"""
Create Snp plots
nknguyen at soe dot ucsc dot edu
Oct 26 2011
"""
import os, sys,re
from optparse import OptionParser
import xml.etree.ElementTree as ET

#from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
import matplotlib.pylab as pylab
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

class Sample():
    def __init__( self, line ):
        items = line.strip().split('\t')
        self.name = items[1]
        self.total = int(items[2])
        self.tp = int(items[5])
        self.fn = 0
        self.sampleTp = self.tp
        
        if len(items) >= 14:
            self.sampleTp = int( items[10] )

    def addsample( self, other ):
        self.total += other.total
        self.tp += other.tp
        self.sampleTp += other.sampleTp
    
    def calcFn(self):
        d = { 'NA19239':[370, 697],
              'NA12878':[417, 739],
              'NA12892':[343, 705],
              'NA19238':[331, 592],
              'NA19240':[277, 595],
              'venter':[319, 270],
              'nigerian':[324, 763],
              'cox':[2393],
              'apd':[683],
              'qbl':[2360],
              'ssto':[2300],
              'dbb':[1975],
              'mann':[1654],
              'mcf':[1545]
            }
        if self.name in d:
            sampleIndels = sum( d[self.name] )
            self.fn = ( sampleIndels - self.sampleTp )*100.0/sampleIndels

    def convertTpToRate(self):
        self.tp = self.tp*100.0/self.total

def readfile( file ):
    f = open(file, 'r')
    samples = {}
    for line in f:
        if re.match('Insertions', line) or re.match('Type', line):
            continue
        items = line.strip().split('\t')
        #if len(items) < 16:
        #    continue
        sample = Sample(line)
        if sample.name not in samples:
            samples[sample.name] = sample
        else:
            samples[sample.name].addsample( sample )
    for s in samples:
        samples[s].calcFn()
        samples[s].convertTpToRate()
    return samples

def readfiles( options ):
    file2expname = { 'indelCheck__0_ignoreAdjacencies_hg19_withAggregates.txt': 'All',\
              'indelCheck__5_ignoreAdjacencies_hg19_withAggregates.txt':'Wobble',\
              'indelCheck_noRepeats_0_ignoreAdjacencies_hg19_withAggregates.txt': 'No-repeats',\
              'indelCheck_noRepeats_5_ignoreAdjacencies_hg19_withAggregates.txt': 'Wobble, No-repeats'}

    exps = {}
    for f in options.files:
        expname = file2expname[ os.path.basename(f) ]
        exps[expname] = readfile( f )
    return exps

def initOptions( parser ):
    parser.add_option( '--outdir', dest='outdir', default='.', help='Output directory' ) 
    parser.add_option( '--indir', dest='indir', help='Required argument. Input directory' ) 
    parser.add_option( '--numOutliners', dest='numOutliners', default=1, help='Number of outliners' ) 
    #parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def checkOptions( args, options, parser ):
    files = [
             'indelCheck__0_ignoreAdjacencies_hg19_withAggregates.txt',\
             'indelCheck__5_ignoreAdjacencies_hg19_withAggregates.txt',\
             'indelCheck_noRepeats_0_ignoreAdjacencies_hg19_withAggregates.txt',\
             'indelCheck_noRepeats_5_ignoreAdjacencies_hg19_withAggregates.txt'
            ]
    if options.indir == None:
        parser.error( 'Please provide input directory\n')
    
    options.files = []
    if not os.path.exists(options.indir):
        parser.error("Input directory %s does not exist\n" %options.indir)

    for f in files:
        filepath = os.path.join( options.indir, f )
        if not os.path.exists( filepath ):
            parser.error( 'File %s does not exist.\n'  % filepath )
        options.files.append( filepath )
    #if options.filteredSamples:
    #    options.filteredSamples = options.filteredSamples.split('-')
    #else:
    #    options.filteredSamples = []

def getSamplesOrder(samples, type):
    nameNstat = []
    rv = True
    for name in samples:
        if type == 'fn' and name == 'yanhuang':
            continue
        if name == 'aggregate' or name == 'panTro3' or name == 'reference':
            continue
        if type == 'tp':
            nameNstat.append( (name, samples[name].tp) )
        elif type == 'fn':
            nameNstat.append( (name, samples[name].fn) )
            rv=False
        elif type == 'total':
            nameNstat.append( (name, samples[name].total) )

    nameNstat = sorted(nameNstat, key=lambda item: item[1], reverse=rv)
    sortedsamples = [item[0] for item in nameNstat]
    #sortedsamples.append('aggregate')
    return sortedsamples

def getData(samples, exps, type):
    data = {} #key = expname, val = list of values cross samples
    for expname in exps:
        data[expname] = []
        currsamples = exps[expname]
        for s in samples:
            if type =='fn' and (s == 'reference' or s =='panTro3' or s=='yanhuang'):
                continue
            if s == 'average':
                data[expname].append( sum(data[expname])/len(data[expname]) )
                continue

            sample = currsamples[s]
            if type == 'tp':
                data[expname].append( sample.tp )
            elif type == 'fn':
                data[expname].append( sample.fn )
            elif type == 'total':
                data[expname].append( sample.total )
    miny = float('inf')
    maxy = float('-inf')
    for exp in data:
        miny = min([ miny, min( data[exp] ) ])
        maxy = max([ maxy, max( data[exp] ) ])
    return data, miny, maxy

def drawPlot(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    #Set title:
    titleDict = {'tp':'Indels confirmed by dbSNP', 'fn':'False Negative', 'total':'Total Indels Called'}
    axes.set_title( titleDict[type] )
    if 'All' not in exps:
        return
    samples = getSamplesOrder( exps['All'], type ) 
    if len( samples ) < 1:
        return

    samples.append('average')
    if type != 'fn':
        samples.append('reference')
        samples.append('panTro3')

    xdata = range( 0, len(samples) )
    #colors = libplot.getColors6()
    colors = libplot.getColors0()
    c = -1
    lines = []
    ydataList, ymin, ymax = getData(samples, exps, type)
    
    scale = -1
    if ymin > 1000:
        scale = len( str(int(ymin)) ) -1
    if scale > 0:
        for exp in ydataList:
            ydataList[exp] = [ float(y)/10**scale for y in ydataList[exp]]
    
    #exporder = ['All', 'Filtered', 'No-repeats', 'Filtered', 'Filtered, No-repeats', 'Recurrent']
    exporder = ['All', 'Wobble', 'No-repeats', 'Wobble, No-repeats']
    if type == 'fn':
        exporder = ['All', 'Wobble']
    elif type == 'total':
        exporder = ['All', 'No-repeats']
    offset = 0.15
    for i, exp in enumerate(exporder):
        xdatai = [x + offset*i for x in xdata]
        ydata = ydataList[exp]
        c += 1
        l = axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=10.0, linestyle='none')
        lines.append(l)

    xmin = -0.4
    xmax = len(samples) - 1 + offset*len(exps) + offset*3
    
    fontP = FontProperties()
    fontP.set_size('x-small')

    if scale > 0:
        ymin = float(ymin)/10**scale
        ymax = float(ymax)/10**scale
    datarange = ymax -ymin
    ymin = ymin - datarange*0.01
    ymax = ymax + datarange*0.01
    
    #Draw vertical lines to separate each sample:
    for i in xrange(1, len(samples)):
        d = (1 - offset*len(exporder))/2.0
        x = [i - d, i - d]
        y = [ymin , ymax]
        axes.plot(x,y, color="#CCCCCC", linewidth=0.005)
    
    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax)
    libplot.editSpine( axes )
 
    axes.set_xticks( [ i + offset*(len(exps)/2-1) for i in range(0, len(samples))] )
    axes.set_xticklabels( [ libplot.properName(s) for s in samples] )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation(90)
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    legend = pyplot.legend( lines, exporder, numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    #axes.set_ylabel( 'Percentage' )
    ylabel = "Percentage"
    if type == 'total':
        ylabel = 'Number of bases'
    if scale > 0:
        ylabel += '(x%d)' %10**scale
    axes.set_ylabel(ylabel)

    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )
   
def drawPlots( options, exps ):
    tpfile = os.path.join( options.outdir, 'indelCheck_TP')
    drawPlot(exps, options, tpfile, "tp")

    fnfile = os.path.join( options.outdir, 'indelCheck_FN' )
    drawPlot(exps, options, fnfile, 'fn')

    totalfile = os.path.join( options.outdir, 'indelCheck_total' )
    drawPlot(exps, options, totalfile, 'total')

def main():
    usage = ('Usage: %prog [options] file1.xml file2.xml\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    exps = readfiles( options )
    drawPlots( options, exps )
    

if __name__ == "__main__":
    main()