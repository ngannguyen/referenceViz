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
                total = int( items[ fields['Total'] ] )
                mapped = int( (items[ fields['Mapped'] -1].split())[0] )
                umapped = int( (items[ fields['UniquelyMapped'] -1 ].split())[0] )
                ppaired = int( (items[ fields['ProperlyPaired'] -1 ].split())[0] )
                uppaired = int( (items[ fields['WithItselftAndMateUniquelyMappedAndProperlyPaired'] -1 ].split())[0] )
                #END HACK

                #total = int( items[ fields['Total'] ] )
                #mapped = int( (items[ fields['Mapped'] ].split())[0] )
                #umapped = int( (items[ fields['UniquelyMapped'] ].split())[0] )
                #ppaired = int( (items[ fields['ProperlyPaired'] ].split())[0] )
                #uppaired = int( (items[ fields['WithItselftAndMateUniquelyMappedAndProperlyPaired'] ].split())[0] )
                
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

def getStats(indir):
    dirs = os.listdir( indir )
    exps = {} #key = sample, val = list of Experiment instances
    for d in dirs:
        dpath = os.path.join( indir, d )
        if not os.path.isdir( dpath ):
            continue
        for r in os.listdir( dpath ):
            rpath = os.path.join(dpath, r)
            if not os.path.isdir( rpath ):
                continue

            exp = Experiment(d, r)

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
    #Get aggregate:
    samples = exps.keys()
    exps['aggregate'] = []
    for s in samples:
        for i, exp in enumerate( exps[s] ):
            aggExpList = exps['aggregate']
            if len(aggExpList) < i+1:
                aggExp = copy.copy(exp)
                #aggExp = copy(exp)
                #aggExp.weight = 1
                aggExp.sample = 'aggregate'
                #aggExp.ref = 'cactusRef'
                exps['aggregate'].append( aggExp )
            else:
                #aggExpList[i].weight += 1
                aggExpList[i].total += exp.total
                aggExpList[i].single += exp.single
                aggExpList[i].paired += exp.paired
                aggExpList[i].mapped += exp.mapped
                aggExpList[i].uniquelyMapped += exp.uniquelyMapped
                aggExpList[i].properlyPaired += exp.properlyPaired
                aggExpList[i].uniquelyMappedAndProperlyPaired += exp.uniquelyMappedAndProperlyPaired
                aggExpList[i].snps += exp.snps
    for exp in exps['aggregate']:
        exp.total /= len(samples)
        exp.single /= len(samples)
        exp.paired /=len(samples)
        exp.mapped /= len(samples)
        exp.uniquelyMapped /= len(samples)
        exp.properlyPaired /= len(samples)
        exp.uniquelyMappedAndProperlyPaired /= len(samples)
        exp.snps /= len(samples)
    #print "AGGREGATE: %d, %s\n" %(len(exps['aggregate']), ",".join( "%s%d" %(e.ref, e.weight) for e in exps['aggregate']) )
                
    return exps

#================= PLOTS ====================
def getData(samples, exps, name):
    #print "getData: samples: %s" %(','.join(samples))
    data = {} #key = refName, val = list of values cross samples
    for sample in samples:
        expList = exps[sample]
        #print "Sample %s, lenExpList %d" %(sample, len(expList))
        for e in expList:
            refname = "%s%d" %(e.ref, e.weight)
            if e.ref == 'hg19':
                refname = 'hg19'
            if refname not in data:
                data[refname] = []
            if name == 'mapped':
                data[refname].append( e.mapped )
            elif name == 'uniquelyMapped':
                data[refname].append( e.uniquelyMapped )
            elif name == 'properlyPaired':
                data[refname].append( e.properlyPaired )
            elif name == 'uniquelyMappedAndProperlyPaired':
                data[refname].append( e.uniquelyMappedAndProperlyPaired )
            elif name == 'snps':
                data[refname].append( e.snps )
    return data

def drawPlot(exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.1, 0.85, 0.85] )

    #Set title:
    titleDict = {'mapped':'Mapped reads', 'uniquelyMapped':'Uniquely Mapped Reads', 'properlyPaired':'Properly Paired Reads', 'uniquelyMappedAndProperlyPaired':'Uniquely Mapped And Properly Paired Reads', 'snps':'Snps'}
    axes.set_title( titleDict[type] )
    
    sampleNhg19mapped = []
    for sample in exps:
        if sample == 'aggregate':
            continue
        explist = exps[sample]
        for e in explist:
            if e.ref == 'hg19':
                sampleNhg19mapped.append( (sample, e.mapped) )

    sampleNhg19mapped = sorted( sampleNhg19mapped, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNhg19mapped]
    samples.append( 'aggregate' )

    xdata = range( 0, len(samples) )
    colors = libplot.getColors0()
    c = -1
    lines = []
    ydataList = getData(samples, exps, type)
    
    refs = sorted( ydataList.keys() )
    miny = float('inf')
    maxy = 0
    offset = 0.05
    axes.set_yscale('log')
    for i, ref in enumerate(refs):
        xdata = [ x + offset for x in xdata ]
        ydata = ydataList[ref]
        miny = min( [miny, min(ydata)] )
        maxy = max( [maxy, max(ydata)] )

        c += 1
        l = axes.plot( xdata, ydata, color=colors[c], marker='.', markersize=10.0, linestyle='none')
        lines.append(l)

    fontP = FontProperties()
    fontP.set_size('x-small')
    axes.set_xlim(-0.5, len(samples) -0.5)
    miny = log10(miny)*0.99
    maxy = log10(maxy)*1.01
    axes.set_ylim( miny, maxy )
    axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.editSpine( axes )

    axes.set_xticks( [ i + offset*(len(refs)/2) for i in range(0, len(samples))] )
    axes.set_xticklabels( samples )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation(90)
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    legend = pyplot.legend( lines, refs, numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    axes.set_ylabel( 'Event counts' )
    libplot.writeImage( fig, pdf, options )

def drawPlots(options, outdir, exps):
    mappedFile = os.path.join(outdir, 'mappingStats-mapped')
    drawPlot(exps, options, mappedFile, "mapped")

    umappedFile = os.path.join(outdir, "mappingStats-uniqMapped.pdf")
    drawPlot(exps, options, umappedFile, "uniquelyMapped")

    ppairedFile = os.path.join(outdir, "mappingStats-properlyPaired.pdf")
    drawPlot(exps, options, ppairedFile, 'properlyPaired')

    uppairedFile = os.path.join(outdir, "mappingStats-uniqProperlyPaired.pdf")
    drawPlot(exps, options, uppairedFile, 'uniquelyMappedAndProperlyPaired')

    snpFile = os.path.join(outdir, "mappingStats-snp.pdf")
    drawPlot(exps, options, snpFile, 'snps')

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
        for e in expList:
            if e.ref == 'hg19':
                f.write("& %s & %s & %s & %s & %s & %s \\\\\n" %(e.ref, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
            else:
                f.write("& %s%d & %s & %s & %s & %s & %s \\\\\n" %(e.ref, e.weight, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
        f.write("\\hline\n")

def makeLatexTab( outfile, exps ):
    f = open(outfile, 'w')
    libplot.writeDocumentStart( f )
    title = "Mapping Stats"
    tabheader ( f, title )
    samples = []
    for s in exps:
        if s != 'aggregate':
            samples.append(s)
    samples.sort()
    samples.append('aggregate')

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

    exps = getStats(options.indir)    
    
    makeLatexTab( os.path.join(options.outdir, "mappingStats.tex"), exps )
    drawPlots(options, options.outdir, exps)


if __name__ == '__main__':
    main()
