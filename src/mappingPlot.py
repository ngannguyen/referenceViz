#!/usr/bin/env python

"""
Summarize mapping stats into tables & plots
Input directory has the structure:
indir/
    experiment/
        catusRef/
            illumina/
                paired/
                    snpcount.txt
                single/
                    snpcount.txt
            summaryStats.txt
        hg19/
            same structure with CactusRef

Output:
"""

import os, sys, re, copy
from optparse import OptionParser
from math import *


#from numpy import *
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.font_manager import FontProperties

from sonLib.bioio import system

class Experiment():
    def __init__( self, expname, ref):
        items = expname.split('_')
        if len(items) < 3:
            sys.stderr.write('Experiment name has unexpected format\n')
            return None
        self.sample = items[2]
        l = items[1].split('-')
        if len(l) < 3:
            sys.stderr.write('Experiment name has unexpected format\n')
            return None
        self.weight = int(l[2] )
        self.ref = ref
        #Initializing the stats:
        self.total = 0
        self.single = 0
        self.paired = 0
        self.mapped = 0
        self.uniquelyMapped = 0
        self.properlyPaired = 0
        self.uniquelyMappedAndProperlyPaired = 0
        self.snps = 0

    def __cmp__(self, other):
        if self.ref < other.ref:
            return -1
        elif self.ref > other.ref:
            return 1
        else:
            if self.weight < other.weight:
                return -1
            elif self.weight == other.weight:
                return 0
            else:
                return 1
    
    def updateStats( self, file ):
        #sys.stderr.write("\nupdateStats, %s, %s, %d, %d, %d, %d, %d. File %s\n" %(self.sample, self.ref, self.weight, self.mapped, self.uniquelyMapped, self.properlyPaired, self.uniquelyMappedAndProperlyPaired, file))
        if not os.path.exists(file):
            return
        f = open(file, 'r')
        fields = {}
        for l in f:
            #sys.stderr.write('%s' %l)
            items = l.strip().split('\t')
            #items = l.strip().split()
            if len(fields) == 0:
                for i, item in enumerate(items):
                    fields[ item.strip() ] = i
            elif len(items) > 0 and items[0] != 'Total':
                #HACK
                #total = int( items[ fields['Total'] ] )
                #mapped = int( (items[ fields['Mapped'] -1].split())[0] )
                #umapped = int( (items[ fields['UniquelyMapped'] -1 ].split())[0] )
                #ppaired = int( (items[ fields['ProperlyPaired'] -1 ].split())[0] )
                #uppaired = int( (items[ fields['WithItselftAndMateUniquelyMappedAndProperlyPaired'] -1 ].split())[0] )
                #END HACK

                total = int( items[ fields['Total'] ] )
                mapped = int( (items[ fields['Mapped'] ].split())[0] )
                umapped = int( (items[ fields['UniquelyMapped'] ].split())[0] )
                ppaired = int( (items[ fields['ProperlyPaired'] ].split())[0] )
                uppaired = int( (items[ fields['WithItselftAndMateUniquelyMappedAndProperlyPaired'] ].split())[0] )
                
                #sys.stderr.write("mapped %d, umapped: %d, ppaired %d, uppaired %d\n" %(mapped, umapped, ppaired, uppaired))
                if re.search('single', items[0]):
                    self.single += total 
                elif re.search('paired', items[0]):
                    self.paired += total
                    self.properlyPaired += ppaired
                    self.uniquelyMappedAndProperlyPaired += uppaired
                self.total += total
                self.mapped += mapped
                self.uniquelyMapped += umapped
        f.close()
        #sys.stderr.write("Done: %s, %s, %d, %d, %d, %d, %d\n" %(self.sample, self.ref, self.weight, self.mapped, self.uniquelyMapped, self.properlyPaired, self.uniquelyMappedAndProperlyPaired))

    def updateSnpCount( self, file ):
        f = open(file, 'r')
        for l in f:
            self.snps += int(l.strip())
        f.close()

def getExpSets(indir):
    dirs = os.listdir( indir )
    patterns = {'required-species-single-copy-species_1-maxWeight-.+-0-True-10000-100-1e-05-False-True-0.0_NA':'True-0.0'
                #'required-species-single-copy-species_1-maxWeight-.+-0-True-10000-100-1e-05-False-False-0.0_NA':'False-0.0',\
                #'required-species-single-copy-species_1-maxWeight-.+-0-True-10000-100-1e-05-False-True-0.2_NA':'True-0.2', \
                #'required-species-single-copy-species_1-maxWeight-.+-0-True-10000-100-1e-05-False-False-0.2_NA':'False-0.2'
                }
    expSets = []
    names = []
    for i, pattern in enumerate( patterns.keys() ):
        expSets.append( [] )
        names.append( patterns[pattern] )
        for d in dirs:
            if re.search( pattern, d ):
                expSets[i].append(d)
    return expSets, names
                
def getStats(indir, dirs):
    #dirs = os.listdir( indir )
    exps = {} #key = sample, val = list of Experiment instances
    for d in dirs:
        if re.search('NA12891', d):#HACK
            continue

        dpath = os.path.join( indir, d )
        if not os.path.isdir( dpath ):
            continue
        for r in os.listdir( dpath ):
            rpath = os.path.join(dpath, r)
            if not os.path.isdir( rpath ):
                continue

            exp = Experiment(d, r)
            if exp.ref == 'cactusRef' and exp.weight == 1: #HACK
                continue

            sumFile = os.path.join( rpath, 'summaryStats.txt' )
            snpFiles = []
            for tech in os.listdir( rpath ):
                techpath = os.path.join( rpath, tech )
                if os.path.isdir( techpath ):
                    for type in os.listdir( techpath ):
                        if type in ['single', 'paired']:
                            snpFiles.append( os.path.join(techpath, type, 'snpcount.txt') )
            
            #if exp.sample in exps and (exp.ref, exp.weight) in [ (e.ref, e.weight) for e in exps[exp.sample] ]:
            if exp.sample in exps and exp.ref in [ e.ref for e in exps[exp.sample] ] and exp.ref == 'hg19':
                continue
            else:
                exp.updateStats( sumFile )
                for snpFile in snpFiles:
                    exp.updateSnpCount( snpFile )
                if exp.sample not in exps:
                    exps[ exp.sample ] = [exp]
                else:
                    exps[ exp.sample ].append( exp )
    #Get average:
    samples = copy.copy( exps.keys() )
    exps['average'] = []
    for s in samples:
        exps[s].sort()
        for i, exp in enumerate( exps[s] ):
            aggExpList = exps['average']
            if len(aggExpList) < i+1:
                aggExp = copy.copy(exp)
                aggExp.sample = 'average'
                exps['average'].append( aggExp )
            else:
                aggExpList[i].total += exp.total
                aggExpList[i].single += exp.single
                aggExpList[i].paired += exp.paired
                aggExpList[i].mapped += exp.mapped
                aggExpList[i].uniquelyMapped += exp.uniquelyMapped
                aggExpList[i].properlyPaired += exp.properlyPaired
                aggExpList[i].uniquelyMappedAndProperlyPaired += exp.uniquelyMappedAndProperlyPaired
                aggExpList[i].snps += exp.snps

    for exp in exps['average']:
        exp.total /= len(samples)
        exp.single /= len(samples)
        exp.paired /=len(samples)
        exp.mapped /= len(samples)
        exp.uniquelyMapped /= len(samples)
        exp.properlyPaired /= len(samples)
        exp.uniquelyMappedAndProperlyPaired /= len(samples)
        exp.snps /= len(samples)
                
    return exps

#================= PLOTS ====================
def getData(samples, exps, name):
    #print "getData: samples: %s" %(','.join(samples))
    data = {} #key = refName, val = list of values cross samples
    
    for sample in samples:
        expList = exps[sample]
        hg19exp = None
        for e in expList:
            if e.ref == 'hg19':
                hg19exp = e
                break

        for e in expList:
            refname = "%s%d" %(e.ref, e.weight)
            if e.ref == 'hg19':
                continue

            if refname not in data:
                data[refname] = []
            
            hg19val = 0
            val = 0
            if name == 'mapped':
                val = e.mapped
                hg19val = hg19exp.mapped
            elif name == 'uniquelyMapped':
                val = e.uniquelyMapped
                hg19val = hg19exp.uniquelyMapped
            elif name == 'properlyPaired':
                val = e.properlyPaired
                hg19val = hg19exp.properlyPaired
            elif name == 'uniquelyMappedAndProperlyPaired':
                val = e.uniquelyMappedAndProperlyPaired
                hg19val = hg19exp.uniquelyMappedAndProperlyPaired
            elif name == 'snps':
                val = e.snps
                hg19val = hg19exp.snps
            
            if hg19val == 0:
                sys.stderr.write('statistics for hg19 is 0.\n')
                data[refname].append( 0 )
            else:
                data[refname].append( (val - hg19val)/float(hg19val) )
    
    miny = float('inf')
    maxy = float('-inf')
    for ref in data:
        currmin = min( data[ref] )
        currmax = max( data[ref] )
        if currmin < miny:
            miny = currmin
        if currmax > maxy:
            maxy = currmax

    return data, miny, maxy

def getYticks(miny, maxy):
    yticks = []
    ylabels = []

    scale = 10**( len(str(miny)) - 1 )
    tick = miny/scale*scale
    yticks.append(tick)
    label = "%.2f" %( float(tick)/scale )
    ylabels.append(label)
    while tick <= maxy:
        tick += scale
        yticks.append(tick)
        label = "%.2f" %( float(tick)/scale )
        ylabels.append(label)

    return yticks, ylabels

def drawSamplePlot(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.14, 0.12, 0.8, 0.8] )

    #Set title:
    axes.set_title( "SNP rate using bwa mapping" )
    
    sampleNsize = []
    for sample in exps:
        if sample == 'average':
            continue
        explist = exps[sample]
        for e in explist:
            if e.ref == 'hg19':
                sampleNsize.append( (sample, e.snps) )

    sampleNsize = sorted( sampleNsize, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNsize]
    samples.append( 'average' )

    #Get ydata:
    ydata1 = [] #hg19
    ydata2 = [] #cactusRef2
    for sample in samples:
        explist = exps[sample]
        for e in explist:
            if e.ref == 'cactusRef' and e.weight == 2:
                ydata2.append( e.snps )
            elif e.ref == 'hg19':
                ydata1.append( e.snps )

    miny = min([min(ydata1), min(ydata2)])
    maxy = max([max(ydata1), max(ydata2)])

    xdata = range( 0, len(samples) )
    colors = ["#E31A1C", "#1F78B4"] #red, blue
    scale = -1
    if miny > 1000:
        scale = len( str(int(miny)) ) - 1
    if scale > 0:
        ydata1 = [ float(y)/10**scale for y in ydata1 ]
        ydata2 = [ float(y)/10**scale for y in ydata2 ]
    lines = []
    lines.append( axes.plot(xdata, ydata1, color=colors[0], marker=".", markersize=16.0, linestyle='none') )
    lines.append( axes.plot(xdata, ydata2, color=colors[1], marker=".", markersize=16.0, linestyle='none') )
    
    if scale > 0:
        miny = float(miny)/10**scale
        maxy = float(maxy)/10**scale

    fontP = FontProperties()
    fontP.set_size('x-small')
    axes.set_xlim(-0.4, len(samples) - 0.6 )
    
    yrange = maxy - miny
    miny = miny - yrange*0.05
    maxy = maxy + yrange*0.1
    axes.set_ylim( miny, maxy )

    libplot.editSpine( axes )

    axes.set_xticks( xdata )
    axes.set_xticklabels( samples )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation(90)
    axes.yaxis.set_ticks_position( 'left' )
    axes.xaxis.set_ticks_position( 'bottom' )
    
    legend = pyplot.legend( lines, [libplot.properName("hg19"), libplot.properName("cactusRef")], numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    axes.set_ylabel( 'Snp counts' )
    if scale > 0:
        axes.set_ylabel( 'Snp counts (x%d)' %(10**scale) )
    axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )

def drawPlot(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    #Set title:
    titleDict = {'mapped':'Mapped reads', 'uniquelyMapped':'Uniquely Mapped Reads', 'properlyPaired':'Properly Paired Reads', 'uniquelyMappedAndProperlyPaired':'Uniquely Mapped And Properly Paired Reads', 'snps':'Snps'}
    axes.set_title( titleDict[type] )
    
    sampleNhg19mapped = []
    for sample in exps:
        if sample == 'average':
            continue
        explist = exps[sample]
        for e in explist:
            if e.ref == 'hg19':
                sampleNhg19mapped.append( (sample, e.total) )

    sampleNhg19mapped = sorted( sampleNhg19mapped, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNhg19mapped]
    samples.append( 'average' )

    xdata = range( 0, len(samples) )
    colors = libplot.getColors0()
    #c = -1
    c = 0
    lines = []
    ydataList, miny, maxy = getData(samples, exps, type)
    
    refs = sorted( ydataList.keys() )
    #miny = float('inf')
    #maxy = 0
    #offset = 0.075
    offset = 0.12
    #axes.set_yscale('log')
    scale = -1
    if miny > 1000:
        scale = len( str(int(miny)) ) - 1
    for i, ref in enumerate(refs):
        xdatai = [ x + offset*i for x in xdata ]
        ydata = ydataList[ref]
        if scale > 0:
            ydata = [ float(y)/10**scale for y in ydata ]
        #miny = min( [miny, min(ydata)] )
        #maxy = max( [maxy, max(ydata)] )
        
        c += 1
        l = axes.plot( xdatai, ydata, color=colors[c], marker='.', markersize=16.0, linestyle='none')
        lines.append(l)

    if scale > 0:
        miny = float(miny)/10**scale
        maxy = float(maxy)/10**scale

    #Draw horizontal line at y = 0:
    xmin = -0.4
    xmax = len(samples) - 1 + offset*len(refs) + offset
    axes.plot( [xmin, xmax], [0,0], color="#6B6B6B", linewidth=0.005)

    fontP = FontProperties()
    fontP.set_size('x-small')
    
    yrange = maxy - miny
    miny = miny - yrange*0.05
    maxy = maxy + yrange*0.2
    #if miny < 0 and maxy >0:
    #    m = max([abs(miny), abs(maxy)])
    #    miny = -m
    #    maxy = m
    
    #Draw vertical lines to separate each sample:
    for i in xrange(1, len(samples)):
        d = (1 - offset*len(refs))/2.0
        x = [i - d, i - d]
        y = [miny , maxy]
        axes.plot(x,y, color="#CCCCCC", linewidth=0.005)
    
    axes.set_xlim(xmin, xmax )
    axes.set_ylim( miny, maxy )
    libplot.editSpine( axes )

    axes.set_xticks( [ i + offset*(len(refs)/2-1) for i in range(0, len(samples))] )
    axes.set_xticklabels( samples )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation(90)
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    properRefs = []
    for r in refs:
        if re.search('cactusRef', r):
            r = r.lstrip('cactusRef')
            properRefs.append( "%s %s" %(libplot.properName('cactusRef'), r))
        else:
            properRefs.append( libplot.properName(r) )

    legend = pyplot.legend( lines,properRefs, numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    axes.set_ylabel( 'Event counts' )
    if scale > 0:
        axes.set_ylabel( 'Event counts (x%d)' %(10**scale) )
    #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )

def drawPlots(options, outdir, exps):
    mappedFile = os.path.join(outdir, 'mappingStats-mapped')
    drawPlot(exps, options, mappedFile, "mapped")

    umappedFile = os.path.join(outdir, "mappingStats-uniqMapped")
    drawPlot(exps, options, umappedFile, "uniquelyMapped")

    ppairedFile = os.path.join(outdir, "mappingStats-properlyPaired")
    drawPlot(exps, options, ppairedFile, 'properlyPaired')

    uppairedFile = os.path.join(outdir, "mappingStats-uniqProperlyPaired")
    drawPlot(exps, options, uppairedFile, 'uniquelyMappedAndProperlyPaired')

    snpFile = os.path.join(outdir, "mappingStats-snp")
    drawPlot(exps, options, snpFile, 'snps')

    file = os.path.join(outdir, "mappingStats-snpSingle")
    drawSamplePlot(exps, options, file, 'snps')

#================= LATEX TABLE ===============
def tabheader(f, title):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{0.75}{%\n")
    f.write("\\begin{tabular}{c|c|r|r|r|r|r}\n")
    f.write("\\multicolumn{7}{c}{%s}\\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("Sample & Ref & Mapped & UniqMapped & PropPaired & UniqPropPaired & Snp\\\\\n")
    f.write("\\hline\n")

def tab( f, exps, samples ):
    for sample in samples:
        expList = exps[sample]
        expList.sort()
        #sys.stderr.write('expList for sample %s: %s\n' %(sample, '\t'.join([ '%s%d' %(e.ref, e.weight)for e in expList])))
        
        f.write( "\\multirow{%d}{*}{%s} " %( len(expList), sample ) )
        #f.write( "\\multirow{%d}{*}{%s} " %( len(expList) -1, sample ) ) #HACK
        for e in expList:
            #if e.ref == 'cactusRef' and e.weight == 1: #HACK
            #    continue
            ref = libplot.properName(e.ref)
            if re.search('cactusRef', e.ref):
                r = e.ref.lstrip('cactusRef')
                ref = "%s %s" % (libplot.properName('cactusRef'), r)

            if e.ref == 'hg19':
                f.write("& %s & %s & %s & %s & %s & %s \\\\\n" %(ref, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
            elif e.ref == 'cactusRef' and e.weight == 2:
                f.write("& \\cellcolor{cyan!30} %s%d & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s \\\\\n" %(ref, e.weight, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
            else:
                f.write("& %s%d & %s & %s & %s & %s & %s \\\\\n" %(ref, e.weight, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
        f.write("\\hline\n")

def makeLatexTab( outfile, exps ):
    f = open(outfile, 'w')
    libplot.writeDocumentStart( f )
    title = "Mapping Stats"
    tabheader ( f, title )
    samples = []
    for s in exps:
        if s != 'average':
            samples.append(s)
    samples.sort()
    samples.append('average')

    tab( f, exps, samples )
    captionStr = "" 
    label = ""
    libplot.tableCloser(f, captionStr, label)
    libplot.writeDocumentEnd( f )
    f.close()

#================= END OF MAKING LATEX TABLE ===============

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory. Required argument.')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory. Default = ./mapStats')
    #parser.add_option()

def checkOptions( args, options, parser ):
    if options.indir == None:
        parser.error('No input directory was given.\n')
    if not os.path.exists( options.indir ):
        parser.error('Input directory does not exist.\n')
    if options.outdir == None:
        options.outdir = os.path.join( os.getcwd(), "mapStats" )
    system("mkdir -p %s" % options.outdir)

def main():
    usage = "mappingPlot.py [options]"
    parser = OptionParser(usage = usage)
    initOptions( parser )
    libplot.initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    expSets, setnames = getExpSets(options.indir)
    for i, dirs in enumerate(expSets):
        name = setnames[i]
        exps = getStats(options.indir, dirs)
        
        outdir = os.path.join( options.outdir, name )
        system("mkdir -p %s" %(outdir))
        makeLatexTab( os.path.join(outdir, "mappingStats.tex"), exps )
        drawPlots(options, outdir, exps)

if __name__ == '__main__':
    main()
