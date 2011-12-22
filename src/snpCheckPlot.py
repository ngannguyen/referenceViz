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
              'dbsnpCheck_noRepeats_hg19_withAggregates.txt': 'No repeats',\
              'dbsnpCheck_noRepeats_filtered_hg19_withAggregates.txt': 'Filtered, No repeats',\
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
        if type == 'tp' or type == 'tpfn' or type == 'tp_all' or type == 'tp_all_recurrent':
            nameNstat.append( (name, samples[name].tp) )
        elif type == 'fn' or type == 'fn_all' or type == 'fn_all_recurrent':
            nameNstat.append( (name, samples[name].fn) )
            rv=False
        elif type == 'total' or type == 'total2' or type == 'total_all':
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
            if (type == 'fn' or type == 'fn_all' or type == 'fn_all_recurrent' ) and (s == 'reference' or s == 'panTro3'):
                continue
            if s == 'average':
                data[expname].append( sum(data[expname])/len(data[expname]) )
                continue

            sample = currsamples[s]
            if type == 'tp' or type == 'tp_all' or type == 'tp_all_recurrent':
                data[expname].append( sample.tp )
            elif type == 'fn' or type == 'fn_all' or type == 'fn_all_recurrent':
                data[expname].append( sample.fn )
            elif type == 'total' or type == 'total2' or type == 'total_all':
                data[expname].append( sample.total )
    miny = float('inf')
    maxy = float('-inf')
    for exp in data:
        miny = min([ miny, min( data[exp] ) ])
        maxy = max([ maxy, max( data[exp] ) ])
    return data, miny, maxy

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

def getOutliers(data):
    #data = {} #key = expname, val = list of values cross samples
    vals = []
    for exp in data:
        for v in data[exp]:
            vals.append(v)
    vals.sort()
    valrange = max(vals) - min(vals)
    cutoffdist = 0.3*valrange
    #print vals
    #print cutoffdist

    prev = vals[0]
    index = len(vals)

    for i, v in enumerate(vals):
        if v - prev >= cutoffdist:
            #print v - prev
            index = i
            break
        prev = v

    if index == len(vals):
        return vals, []
    #print vals[:index]
    #print vals[index:]
    return vals[:index], vals[index:]

def drawPlot2(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 11.2, 10.0, options )

    #Set title:
    titleDict = {'total':'Total SNPs Called', 'total2':'Total SNPs Called', 'total_all':'Total SNPs Called'}
    if 'All' not in exps:
        return
    samples = getSamplesOrder( exps['All'], type ) 
    if len( samples ) < 1:
        return
    
    samples.append('average')
    samples.append('reference')
    samples.append('panTro3')
    
    xdata = range( 0, len(samples) )
    colors = libplot.getColors6()
    c = -1
    lines = []
    
    pointsize = 10.0
    offset = 0.15
    exporder = ['All', 'Recurrent', 'No repeats', 'Filtered', 'Filtered, No repeats']
    if type == 'total_all':
        exporder = ['All']
    elif type == 'total2':
        exporder = ['All', 'Recurrent']
     
    #Get ydata
    ydataList, ymin, ymax = getData(samples, exps, type, exporder)
    yrange = ymax - ymin

    #Get normal range and outlier range:
    normalvals, outliers = getOutliers(ydataList)
    minNormal = min(normalvals) - 0.02*yrange
    maxNormal = max(normalvals) + 0.02*yrange 
    minOutlier = min(outliers) - 0.02*yrange
    maxOutlier = max(outliers) + 0.02*yrange
    if minNormal < 0:
        minNormal = -0.5
    
    #Set up the axes
    ax, ax2 = setAxes( fig, maxOutlier - minOutlier, maxNormal - minNormal)
     
    scale = -1
    if minNormal > 5000:
        scale = len( str(int(minNormal)) ) -1
    if scale > 0:
        for exp in ydataList:
            ydataList[exp] = [ float(y)/10**scale for y in ydataList[exp]]
    
    #PLOT 
    for i, exp in enumerate(exporder):
        xdatai = [x + offset*i for x in xdata]
        ydata = ydataList[exp]
        c += 1
        #Outlier plot
        l = ax.plot(xdatai, ydata, color=colors[c], marker='.', markersize=pointsize, linestyle='none')
        #Normal range plot
        ax2.plot(xdatai, ydata, color=colors[c], marker='.', markersize=pointsize, linestyle='none')
        lines.append(l)

    xmin = -0.4
    xmax = len(samples) - 1 + offset*len(exps) + offset*3
    
    fontP = FontProperties()
    fontP.set_size('x-small')

    if scale > 0:
        minNormal = float(minNormal)/10**scale
        maxNormal = float(maxNormal)/10**scale
        minOutlier = float(minOutlier)/10**scale
        maxOutlier = float(maxOutlier)/10**scale
        #ymin = float(ymin)/10**scale
        #ymax = float(ymax)/10**scale
    #datarange = ymax - ymin
    #ymin = ymin - datarange*0.01
    #ymax = ymax + datarange*0.01
    
    #Draw the Discontinue sign:
    d = 0.2 #how big to make the diagonal lines in axes coordinates
    if scale == -1:
        d = 150
    ax.plot( (-1, 0), (minOutlier +d, minOutlier - d), color = "k", clip_on=False )
    ax2.plot( (-1, 0), (maxNormal +d, maxNormal - d), color = "k", clip_on=False )

    #Draw vertical lines to separate each sample:
    for i in xrange(1, len(samples)):
        d = (1 - offset*len(exporder))/2.0
        x = [i - d, i - d]
        y = [minNormal , maxOutlier]
        ax.plot(x,y, color="#CCCCCC", linewidth=0.005)
        ax2.plot(x,y, color="#CCCCCC", linewidth=0.005)
    
    #axes.set_ylim(ymin, ymax)
    #libplot.editSpine( axes )
    xticklabels = [libplot.properName(s) for s in samples]

    #Set limits:
    ax.set_ylim( minOutlier, maxOutlier )
    ax.set_xlim(xmin, xmax)
    ax.set_xticks( [ i + offset*(len(exps)/2-1) for i in range(0, len(samples))] )
    dummyxticklabels = [ "" for l in xticklabels ]
    ax.set_xticklabels(dummyxticklabels)
    
    #Make sure the y ticks of the top plot is the same with the bottom plot:
    step = 2
    if scale == -1:
        step = 5000
    ytickpositions = []
    ytickpos = 0
    while ytickpos < maxOutlier:
        if ytickpos >= minOutlier:
            ytickpositions.append(ytickpos)
        ytickpos += step
    ax.set_yticks(ytickpositions)

    #Set limit for the bottom plot:
    ax2.set_ylim(minNormal, maxNormal)
    ax2.set_xlim(xmin, xmax)

    #Hide the spines between ax and ax2:
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('none')

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xaxis.tick_bottom()
    ax2.yaxis.set_ticks_position( 'left' )

    ax2.set_xticks( [ i + offset*(len(exps)/2-1) for i in range(0, len(samples))] )
    ax2.set_xticklabels( xticklabels ) 
    
    for label in ax2.xaxis.get_ticklabels():
        label.set_rotation(75)
   
    if len(exporder) > 1:
        legend = pyplot.legend( lines, exporder, numpoints=1, loc='upper left', prop=fontP)
        legend._drawFrame = False

    ax2.set_xlabel( 'Samples' )
    ylabel = 'Number of SNPs'
    if scale > 0:
        ylabel += ' (x%d)' %10**scale
    ax2.set_ylabel(ylabel)
    ax.set_title( titleDict[type] )

    ax.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    ax2.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )

def drawPlot(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 11.2, 10.0, options )
    axes = fig.add_axes( [0.12, 0.18, 0.85, 0.75] )

    #Set title:
    titleDict = {'tp':'True Positives According to dbSNP', 'fn':'False Negatives According to dbSNP', \
                 'tp_all':'True Positives According to dbSNP', 'fn_all':'False Negatives According to dbSNP', \
                 'tp_all_recurrent':'True Positives According to dbSNP', 'fn_all_recurrent':'False Negatives According to dbSNP', \
                 'total':'Total SNPs Called', 'total2':'Total SNPs Called', 'total_all':'Total SNPs Called', \
                 'tpfn':'SNP Overlap with dbSNP'}

    axes.set_title( titleDict[type] )
    if 'All' not in exps:
        return
    samples = getSamplesOrder( exps['All'], type ) 
    if len( samples ) < 1:
        return
    
    samples.append('average')
    if type != 'fn' and type !='fn_all' and type != 'fn_all_recurrent':
        samples.append('reference')
        samples.append('panTro3')

    xdata = range( 0, len(samples) )
    colors = libplot.getColors6()
    c = -1
    lines = []
    
    pointsize = 10.0
    offset = 0.15
    #exporder = ['All', 'Filtered', 'No-repeats', 'Filtered, No-repeats', 'Recurrent']
    exporder = ['All', 'Recurrent', 'No repeats', 'Filtered', 'Filtered, No repeats']
    if type == 'total_all' or type == 'tp_all' or type == 'fn_all':
        exporder = ['All']
        #pointsize = 16.0
        #offset = 0.3
    elif type == 'total2' or type == 'tpfn' or type == 'fn_all_recurrent' or type == 'tp_all_recurrent':
        exporder = ['All', 'Recurrent']
        #pointsize = 16.0
        #offset = 0.3
    
    ydataList, ymin, ymax = getData(samples, exps, type, exporder)
     
    scale = -1
    if ymin > 1000:
        scale = len( str(int(ymin)) ) -1
    if scale > 0:
        for exp in ydataList:
            ydataList[exp] = [ float(y)/10**scale for y in ydataList[exp]]
    
    if type == 'tpfn':
        for j, t in enumerate(['tp', 'fn']):
            for i, exp in enumerate(exporder):
                if t == 'tp':
                    xdatai = [x + offset*(j*2+i) for x in xdata]
                else:
                    xdatai = [x + offset*(j*2+i) for x in xdata[: len(xdata) -2] ]
                ydata = ydataList["%s.%s" %(exp,t)]
                c += 1
                lines.append(axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=pointsize, linestyle='none'))
    else:
        for i, exp in enumerate(exporder):
            xdatai = [x + offset*i for x in xdata]
            ydata = ydataList[exp]
            c += 1
            l = axes.plot(xdatai, ydata, color=colors[c], marker='.', markersize=pointsize, linestyle='none')
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
        legend = pyplot.legend( lines, ['All, TP','Recurrent, TP','All, FN','Recurrent, FN'], numpoints=1, loc='best', prop=fontP)
        legend._drawFrame = False
    elif type == 'tp_all_recurrent' or type == 'tp' or type == 'tp_all':
        legend = pyplot.legend( lines, exporder, numpoints=1, loc='lower left', prop=fontP)
        legend._drawFrame = False
    elif len(exporder) > 1:
        legend = pyplot.legend( lines, exporder, numpoints=1, loc='best', prop=fontP)
        legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    #axes.set_ylabel( 'Percentage' )
    ylabel = "Percentage"
    if type == 'total' or type == 'total2':
        ylabel = 'Number of SNPs'
    if scale > 0:
        ylabel += ' (x%d)' %10**scale
    axes.set_ylabel(ylabel)

    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )
   
def drawPlots( options, exps ):
    tpfile = os.path.join( options.outdir, 'dbsnpCheck_TP')
    drawPlot(exps, options, tpfile, "tp")

    fnfile = os.path.join( options.outdir, 'dbsnpCheck_FN' )
    drawPlot(exps, options, fnfile, 'fn')

    totalfile = os.path.join( options.outdir, 'dbsnpCheck_total' )
    drawPlot2(exps, options, totalfile, 'total')

    total2file = os.path.join( options.outdir, 'dbsnpCheck_total2' )
    drawPlot2(exps, options, total2file, 'total2')

    tpfnfile = os.path.join( options.outdir, 'dbsnpCheck_TPFN' )
    drawPlot(exps, options, tpfnfile, 'tpfn')

    total2file = os.path.join( options.outdir, 'dbsnpCheck_total_all' )
    drawPlot2(exps, options, total2file, 'total_all')
    
    tpallfile = os.path.join( options.outdir, 'dbsnpCheck_TP_all' )
    drawPlot(exps, options, tpallfile, 'tp_all')
    
    fnallfile = os.path.join( options.outdir, 'dbsnpCheck_FN_all' )
    drawPlot(exps, options, fnallfile, 'fn_all')

    tpfnfile = os.path.join( options.outdir, 'dbsnpCheck_TP_all_recurrent' )
    drawPlot(exps, options, tpfnfile, 'tp_all_recurrent')
    
    tpfnfile = os.path.join( options.outdir, 'dbsnpCheck_FN_all_recurrent' )
    drawPlot(exps, options, tpfnfile, 'fn_all_recurrent')

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
