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
        #if len(items) < 15:
            #print line
            #print len(items)
        self.name = items[0]
        self.total = int(items[1])
        self.tp = float(items[5])
        self.fn = 0
        if self.name != 'reference' and self.name != 'panTro3':
            self.fn = float(items[14])
            sample2snps = { #HACK
                            'NA19239' : 15026,
                            'NA12878' : 15108,
                            'NA12892' : 13123,
                            'venter'  : 14523,
                            'yanhuang': 15902,
                            'NA19238' : 13312,
                            'nigerian': 15658,
                            'NA19240' : 14358
                        }
            if self.name in sample2snps:
                self.fn = (sample2snps[self.name] - int(items[11]))*100.0/sample2snps[self.name]

        #self.totalbases = int( items[len(items) -1] )

def readfile( file ):
    f = open(file, 'r')
    samples = {}
    for line in f:
        if re.match('TotalSnps', line) or re.match('Sample', line):
            continue
        items = line.strip().split('\t')
        if not re.match('reference', line) and not re.match('panTro3', line) and len(items) < 16:
            #print "Less than 16 fields: %s" %line
            continue
        sample = Sample(line)
        samples[sample.name] = sample
    
    #Add aggregate:
    #for s in samples:
    #    if s == 'reference' or s == 'panTro3':
    #        continue
    #    sample = samples[s]

    return samples

def readfiles( options ):
    file2expname = { 'dbsnpCheck_hg19_withAggregates.txt': 'All',\
              'dbsnpCheck_filtered_hg19_withAggregates.txt':'Filtered',\
              'dbsnpCheck_noRepeats_hg19_withAggregates.txt': 'No-repeats',\
              'dbsnpCheck_noRepeats_filtered_hg19_withAggregates.txt': 'Filtered, No-repeats',\
              'dbsnpCheck_hg19_recurrent_withAggregates.txt': 'Recurrent'}

    exps = {}
    for f in options.files:
        expname = file2expname[ os.path.basename(f) ]
        exps[expname] = readfile( f )
    return exps

def initOptions( parser ):
    parser.add_option( '--outdir', dest='outdir', default='.', help='Output directory' ) 
    parser.add_option( '--indir', dest='indir', help='Required argument. Input directory' ) 
    parser.add_option( '--numOutliners', dest='numOutliners', default=1, help='Number of outliners' ) 
    parser.add_option('--filteredSamples', dest='filteredSamples', help='Hyphen separated list of samples that were filtered out (not to include in the plot)')

def checkOptions( args, options, parser ):
    files = [ 'dbsnpCheck_hg19_withAggregates.txt',\
              'dbsnpCheck_filtered_hg19_withAggregates.txt',\
              'dbsnpCheck_noRepeats_hg19_withAggregates.txt',\
              'dbsnpCheck_noRepeats_filtered_hg19_withAggregates.txt',\
              'dbsnpCheck_hg19_recurrent_withAggregates.txt']

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
        if name == 'aggregate' or name == 'panTro3' or name == 'reference':
            continue
        if type == 'tp' or type == 'tpfn':
            nameNstat.append( (name, samples[name].tp) )
        elif type == 'fn':
            nameNstat.append( (name, samples[name].fn) )
            rv=False
        elif type == 'total' or type == 'total2':
            nameNstat.append( (name, samples[name].total) )

    nameNstat = sorted(nameNstat, key=lambda item: item[1], reverse=rv)
    sortedsamples = [item[0] for item in nameNstat]
    #sortedsamples.append('average')
    #sortedsamples.append('reference')
    #sortedsamples.append('panTro3')
    return sortedsamples

def getData(samples, exps, type, expnames):
    data = {} #key = expname, val = list of values cross samples
    if type == 'tpfn':
        tpdata,tpmin, tpmax = getData(samples, exps, 'tp', expnames)
        fndata,fnmin, fnmax = getData(samples, exps, 'fn', expnames)
        for exp in tpdata:
            data["%s.tp" %exp] = tpdata[exp]
        for exp in fndata:
            data["%s.fn" %exp] = fndata[exp]
        return data, min([tpmin, fnmin]), max([tpmax, fnmax])
        
    for expname in exps:
        if expname not in expnames:
            continue
        data[expname] = []
        currsamples = exps[expname]
        for s in samples:
            if type == 'fn' and (s == 'reference' or s == 'panTro3'):
                continue
            if s == 'average':
                data[expname].append( sum(data[expname])/len(data[expname]) )
                continue

            sample = currsamples[s]
            if type == 'tp':
                data[expname].append( sample.tp )
            elif type == 'fn':
                data[expname].append( sample.fn )
            elif type == 'total' or type == 'total2':
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
    axes = fig.add_axes( [0.12, 0.18, 0.85, 0.75] )

    #Set title:
    titleDict = {'tp':'Snps confirmed by dbSNP', 'fn':'False Negative', 'total':'Total Snps Called', 'total2':'Total Snps Called', 'tpfn':''}
    axes.set_title( titleDict[type] )
    if 'All' not in exps:
        return
    samples = getSamplesOrder( exps['All'], type ) 
    print samples
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
    
    #exporder = ['All', 'Filtered', 'No-repeats', 'Filtered, No-repeats', 'Recurrent']
    exporder = ['All', 'Recurrent', 'No-repeats', 'Filtered', 'Filtered, No-repeats']
    if type == 'total2' or type == 'tpfn':
        exporder = ['All', 'Recurrent']
    
    ydataList, ymin, ymax = getData(samples, exps, type, exporder)
     
    scale = -1
    if ymin > 1000:
        scale = len( str(int(ymin)) ) -1
    if scale > 0:
        for exp in ydataList:
            ydataList[exp] = [ float(y)/10**scale for y in ydataList[exp]]
    
    offset = 0.15
    for i, exp in enumerate(exporder):
        if type == 'tpfn':
            for j, t in enumerate(['tp', 'fn']):
                if t == 'tp':
                    xdatai = [x + offset*(i*2+j) for x in xdata]
                else:
                    xdatai = [x + offset*(i*2+j) for x in xdata[: len(xdata) -2] ]
                ydata = ydataList["%s.%s" %(exp,t)]
                c += 1
                lines.append(axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=10.0, linestyle='none'))
        else:
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
    datarange = ymax - ymin
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
        label.set_rotation(75)
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
   
    if type == 'tpfn':
        legend = pyplot.legend( lines, ['All, TP','All, FN','Recurrent, TP','Recurrent, FN'], numpoints=1, loc='best', prop=fontP)
    else:
        legend = pyplot.legend( lines, exporder, numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    #axes.set_ylabel( 'Percentage' )
    ylabel = "Percentage"
    if type == 'total':
        ylabel = 'Number of snps'
    if scale > 0:
        ylabel += '(x%d)' %10**scale
    axes.set_ylabel(ylabel)

    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )
   
def drawPlots( options, exps ):
    tpfile = os.path.join( options.outdir, 'dbsnpCheck_TP')
    drawPlot(exps, options, tpfile, "tp")

    fnfile = os.path.join( options.outdir, 'dbsnpCheck_FN' )
    drawPlot(exps, options, fnfile, 'fn')

    totalfile = os.path.join( options.outdir, 'dbsnpCheck_total' )
    drawPlot(exps, options, totalfile, 'total')

    total2file = os.path.join( options.outdir, 'dbsnpCheck_total2' )
    drawPlot(exps, options, total2file, 'total2')

    tpfnfile = os.path.join( options.outdir, 'dbsnpCheck_TPFN' )
    drawPlot(exps, options, tpfnfile, 'tpfn')

def getFilteredSamples(exp):
    #filtered sample:
    filteredSamples = ''
    items = exp.split('_')
    if len(items) > 2:
        filteredSamples = items[ len(items) -1 ]
    return filteredSamples

def main():
    usage = ('Usage: %prog [options] file1.xml file2.xml\n\n')
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )

    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    #HACK:
    filteredSamples = getFilteredSamples( os.path.basename(options.indir) )
    if len(filteredSamples) > 0:
        return

    exps = readfiles( options )
    drawPlots( options, exps )
    

if __name__ == "__main__":
    main()
