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
                merge.bam
                snpcount.txt
                mergeSorted.bam (use this to generate the pileup file, which gives the total number of bases the mapped reads cover)
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
        self.weight = 0
        if len(l) >= 3:
            self.weight = int(l[2])
            #sys.stderr.write('Experiment name has unexpected format\n')
            #return None
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
        self.totalbases = 0

        #Single & paired stats:
        self.singleMapped = 0
        self.singleUniquelyMapped = 0
        self.pairedMapped = 0
        self.pairedUniquelyMapped = 0
        self.pairedProperlyPaired = 0
        self.pairedUniquelyMappedAndProperlyPaired = 0

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
                total = int( items[ fields['Total'] ] )
                mapped = int( (items[ fields['Mapped'] ].split())[0] )
                umapped = int( (items[ fields['UniquelyMapped'] ].split())[0] )
                ppaired = int( (items[ fields['ProperlyPaired'] ].split())[0] )
                uppaired = int( (items[ fields['WithItselftAndMateUniquelyMappedAndProperlyPaired'] ].split())[0] )
                
                if re.search('single', items[0]):
                    self.single += total 
                    self.singleMapped += mapped
                    self.singleUniquelyMapped += umapped

                elif re.search('paired', items[0]):
                    self.paired += total
                    self.properlyPaired += ppaired
                    self.uniquelyMappedAndProperlyPaired += uppaired
                    
                    self.pairedMapped += mapped
                    self.pairedUniquelyMapped += umapped
                    self.pairedProperlyPaired += ppaired
                    self.pairedUniquelyMappedAndProperlyPaired += uppaired
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

    def getTotalBases(self, dir, file):
        pileupFile = os.path.join(dir, "mergePileup.txt")
        totalCountFile = os.path.join(dir, "tempPileupCount.txt")
        if not os.path.exists(pileupFile):
            system("samtools mpileup %s > %s" %(file, pileupFile))

        if not os.path.exists(totalCountFile):
            system("cat %s | grep -vc '^#' > %s" %(pileupFile, totalCountFile))
        f = open(totalCountFile, 'r')
        line = f.readline()
        f.close()
        self.totalbases = int(line.strip())
        self.snprate = 0.0
        if self.totalbases > 0:
            self.snprate = float(self.snps)/self.totalbases
        #sys.stderr.write("%s\t%d\t%f\n" %(dir, self.totalbases, self.snprate))
        #system("rm %s" %totalCountFile)

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
                
def getRefStats(indir):
    exps = {} #key = ref(e.g hg19, apd...), val = {} in which key = sample, val = list of Experiment instances
    samples = os.listdir(indir)
    for sample in samples:
        sampledir = os.path.join(indir, sample)
        if not os.path.isdir( sampledir ):
            continue
        refs = os.listdir( sampledir )
        for ref in refs:
            rpath = os.path.join(sampledir, ref)
            if not os.path.isdir( rpath ):
                continue

            dummyExpname = 'ref_%s_%s' %(ref, sample)
            exp = Experiment(dummyExpname, ref)

            sumFile = os.path.join( rpath, 'summaryStats.txt' )
            snpFiles = []
            bamfile = ''
            tpath = ''
            for tech in os.listdir( rpath ):
                techpath = os.path.join( rpath, tech )
                if os.path.isdir( techpath ):
                    snpFiles.append( os.path.join(techpath, 'snpcount.txt') )
                    bamfile =  os.path.join(techpath, 'mergeSorted.bam')
                    tpath = techpath
                    #for type in os.listdir( techpath ):
                    #    if type in ['single', 'paired']:
                    #        snpFiles.append( os.path.join(techpath, type, 'snpcount.txt') )
                
            exp.updateStats( sumFile )
            for snpFile in snpFiles:
                exp.updateSnpCount( snpFile )
            exp.getTotalBases(tpath, bamfile)
           
            if ref not in exps:
                exps[ref] = { sample:exp }
            else:
                if sample in exps[ref]:
                    sys.stderr.write("Ref %s has repeated sample %s\n" %(ref, sample))
                exps[ref][sample] = exp
    
    #Get average:
    for r in exps.keys():
        rexps = exps[r]
        samples = copy.copy( rexps.keys() )
        for s in samples:
            exp = rexps[s]
            if 'average' not in rexps:
                avrExp = copy.copy( exp )
                avrExp.sample = 'average'
                rexps['average'] = avrExp
            else:
                rexps['average'].total += exp.total
                rexps['average'].single += exp.single
                rexps['average'].paired += exp.paired
                rexps['average'].mapped += exp.mapped
                rexps['average'].uniquelyMapped += exp.uniquelyMapped
                rexps['average'].properlyPaired += exp.properlyPaired
                rexps['average'].uniquelyMappedAndProperlyPaired += exp.uniquelyMappedAndProperlyPaired
                #rexps['average'].snps += exp.snps
                #rexps['average'].totalbases += exp.totalbases
                rexps['average'].snprate += exp.snprate
                
                rexps['average'].singleMapped += exp.singleMapped
                rexps['average'].singleUniquelyMapped += exp.singleUniquelyMapped
                rexps['average'].pairedMapped += exp.pairedMapped
                rexps['average'].pairedUniquelyMapped += exp.pairedUniquelyMapped
                rexps['average'].pairedProperlyPaired += exp.pairedProperlyPaired
                rexps['average'].pairedUniquelyMappedAndProperlyPaired += exp.pairedUniquelyMappedAndProperlyPaired
                        
        rexps['average'].total /= len(samples)
        rexps['average'].single /= len(samples)
        rexps['average'].paired /=len(samples)
        rexps['average'].mapped /= len(samples)
        rexps['average'].uniquelyMapped /= len(samples)
        rexps['average'].properlyPaired /= len(samples)
        rexps['average'].uniquelyMappedAndProperlyPaired /= len(samples)
        #rexps['average'].snps /= len(samples)
        rexps['average'].snprate /= len(samples)
        rexps['average'].singleMapped /= len(samples)
        rexps['average'].singleUniquelyMapped /= len(samples)
        rexps['average'].pairedMapped /= len(samples)
        rexps['average'].pairedUniquelyMapped /= len(samples)
        rexps['average'].pairedProperlyPaired /= len(samples)
        rexps['average'].pairedUniquelyMappedAndProperlyPaired /= len(samples)
    
        sys.stderr.write("%s\n" %r)
        sys.stderr.write('SingleReads: %d\t' %rexps['average'].single)
        sys.stderr.write('singleMapped: %f\t' %rexps['average'].singleMapped)
        sys.stderr.write('uiniquelyMapped: %f\n' %rexps['average'].uniquelyMapped)
        sys.stderr.write('PairedReads: %d\t' %rexps['average'].paired)
        sys.stderr.write('pairedMapped: %f\t' %rexps['average'].pairedMapped)
        sys.stderr.write('pairedUniquelyMapped: %f\t' %rexps['average'].pairedUniquelyMapped)
        sys.stderr.write('pairedProperlyPaired: %f\t' %rexps['average'].pairedProperlyPaired)
        sys.stderr.write('pairedUniquelyProperlyPaired: %f\n' %rexps['average'].pairedUniquelyMappedAndProperlyPaired)

    return exps

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
            #if exp.ref == 'cactusRef' and exp.weight == 1: #HACK
            if exp.ref == 'cactusRef' and exp.weight == 8: #HACK
                continue

            sumFile = os.path.join( rpath, 'summaryStats.txt' )
            snpFiles = []
            bamfile = ''
            tpath = ''
            for tech in os.listdir( rpath ):
                techpath = os.path.join( rpath, tech )
                if os.path.isdir( techpath ):
                    snpFiles.append( os.path.join(techpath, 'snpcount.txt') )
                    bamfile =  os.path.join(techpath, 'mergeSorted.bam')
                    tpath = techpath
                    #for type in os.listdir( techpath ):
                    #    if type in ['single', 'paired']:
                    #        snpFiles.append( os.path.join(techpath, type, 'snpcount.txt') )
            
            #if exp.sample in exps and (exp.ref, exp.weight) in [ (e.ref, e.weight) for e in exps[exp.sample] ]:
            if exp.sample in exps and exp.ref in [ e.ref for e in exps[exp.sample] ] and exp.ref == 'hg19':
                continue
            else:
                exp.updateStats( sumFile )
                for snpFile in snpFiles:
                    exp.updateSnpCount( snpFile )
                exp.getTotalBases(tpath, bamfile)
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
                #aggExpList[i].snps += exp.snps
                aggExpList[i].snprate += exp.snprate
                
                aggExpList[i].singleMapped += exp.singleMapped
                aggExpList[i].singleUniquelyMapped += exp.singleUniquelyMapped
                aggExpList[i].pairedMapped += exp.pairedMapped
                aggExpList[i].pairedUniquelyMapped += exp.pairedUniquelyMapped
                aggExpList[i].pairedProperlyPaired += exp.pairedProperlyPaired
                aggExpList[i].pairedUniquelyMappedAndProperlyPaired += exp.pairedUniquelyMappedAndProperlyPaired
                

    for exp in exps['average']:
        exp.total /= len(samples)
        exp.single /= len(samples)
        exp.paired /=len(samples)
        exp.mapped /= len(samples)
        exp.uniquelyMapped /= len(samples)
        exp.properlyPaired /= len(samples)
        exp.uniquelyMappedAndProperlyPaired /= len(samples)
        #exp.snps /= len(samples)
        exp.snprate /= len(samples)

        exp.singleMapped /= len(samples)
        exp.singleUniquelyMapped /= len(samples)
        exp.pairedMapped /= len(samples) 
        exp.pairedUniquelyMapped /= len(samples)
        exp.pairedProperlyPaired /= len(samples) 
        exp.pairedUniquelyMappedAndProperlyPaired /= len(samples)

        sys.stderr.write("%s, %d\n" %(exp.ref, exp.weight))
        sys.stderr.write('SingleReads: %d\t' %exp.single)
        sys.stderr.write('singleMapped: %f\t' % exp.singleMapped)
        sys.stderr.write('uiniquelyMapped: %f\n' %exp.uniquelyMapped)
        sys.stderr.write('PairedReads: %d\t' %exp.paired)
        sys.stderr.write('pairedMapped: %f\t' %exp.pairedMapped)
        sys.stderr.write('pairedUniquelyMapped: %f\t' %exp.pairedUniquelyMapped)
        sys.stderr.write('pairedProperlyPaired: %f\t' %exp.pairedProperlyPaired)
        sys.stderr.write('pairedUniquelyProperlyPaired: %f\n' %exp.pairedUniquelyMappedAndProperlyPaired)
                
    return exps

#================= PLOTS ====================
def getData(samples, exps, rexps, name):
    #print "getData: samples: %s" %(','.join(samples))
    data = {} #key = refName, val = list of values cross samples
    
    for sample in samples:
        expList = exps[sample]
        hg19exp = rexps[sample]

        for e in expList:
            refname = "%s%d" %(e.ref, e.weight)
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
                #val = e.snps
                #hg19val = hg19exp.snps
                val = e.snprate
                hg19val = hg19exp.snprate
            
            if hg19val == 0:
                sys.stderr.write('%s: statistics for hg19 is 0.\n' %name)
                data[refname].append( 0 )
            else:
                data[refname].append( (val - hg19val)*100.0/float(hg19val) )
    
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

def getData2(samples, rexps, exps, names):
    data = {} #key = , val = list of values cross samples
    
    for name in names:
        data[name] = []
        for sample in samples:
            expList = exps[sample]
            hg19exp = rexps[sample]

            for e in expList:
                refname = "%s%d" %(e.ref, e.weight)
                if refname != 'cactusRef2':
                    continue
                
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
                    #val = e.snps
                    #hg19val = hg19exp.snps
                    val = e.snprate
                    hg19val = hg19exp.snprate
                
                if hg19val == 0:
                    sys.stderr.write('statistics for hg19 is 0.\n')
                    data[name].append( 0 )
                else:
                    data[name].append( (val - hg19val)*100.0/float(hg19val) )
                break
        
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

def drawSamplePlot(rexps, exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 11.2, 10.0, options )
    axes = fig.add_axes( [0.14, 0.12, 0.8, 0.8] )

    #Set title:
    axes.set_title( "SNP Rate Using BWA Mapping" )
    
    sampleNsize = []
    if len(rexps) < 1:
        return
    ref = ''
        
    for sample in rexps:
        if sample == 'average':
            continue
        exp = rexps[sample]
        ref = exp.ref
        #sampleNsize.append( (sample, exp.snps) )
        sampleNsize.append( (sample, exp.snprate) )
    otherRefName = ref

    sampleNsize = sorted( sampleNsize, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNsize]
    samples.append( 'average' )

    #Get ydata:
    ydata1 = [] #otherRef (hg19, apd, ...)
    ydata2 = [] #cactusRef2
    for sample in samples:
        explist = exps[sample]
        otherRef = rexps[sample]
        ydata1.append( otherRef.snprate )
        for e in explist:
            if e.ref == 'cactusRef' and e.weight == 2:
                ydata2.append( e.snprate )

    miny = min([min(ydata1), min(ydata2)])
    maxy = max([max(ydata1), max(ydata2)])

    xdata = range( 0, len(samples) )
    #colors = ["#E31A1C", "#1F78B4"] #red, blue
    colors = ["#1F78B4", "#E31A1C"] #red, blue
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
    
    legend = pyplot.legend( lines, [libplot.properName(otherRefName), libplot.properName("cactusRef")], numpoints=1, loc='best', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    axes.set_ylabel( 'SNPs Per Site' )
    if scale > 0:
        axes.set_ylabel( 'Snp counts (x%d)' %(10**scale) )
    axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )

def drawPlot(rexps, exps, options, outfile, type):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.14, 0.85, 0.8] )

    #Set title:
    titleDict = {'mapped':'Mapped reads', 'uniquelyMapped':'Uniquely Mapped Reads', 'properlyPaired':'Properly Paired Reads', 'uniquelyMappedAndProperlyPaired':'Uniquely Mapped And Properly Paired Reads', 'snps':'SNPs'}
    axes.set_title( titleDict[type] )
    
    if len(rexps) < 1:
        return
    
    sampleNotherRefmapped = []
    ref = ''
    for sample in rexps:
        if sample == 'average':
            continue
        exp = rexps[sample]
        ref = exp.ref
        sampleNotherRefmapped.append( (sample, exp.total) )
    otherRefName = libplot.properName( ref )
    
    sampleNotherRefmapped = sorted( sampleNotherRefmapped, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNotherRefmapped]
    samples.append( 'average' )

    xdata = range( 0, len(samples) )
    colors = libplot.getColors4()
    #c = -1
    c = 0
    lines = []
    ydataList, miny, maxy = getData(samples, exps, rexps, type)
    #print ydataList
    
    refs = sorted( ydataList.keys() )
    #miny = float('inf')
    #maxy = 0
    #offset = 0.075
    offset = 0.12
    #if type != 'snps':
    #    offset = 0
    #axes.set_yscale('log')
    scale = -1
    if miny > 1000:
        scale = len( str(int(miny)) ) - 1

    #Draw line connecting the data for each sample (each bin):
    binXdataList = [ [] for x in xdata ]
    binYdataList = [ [] for x in xdata ]
    for i, ref in enumerate(refs):
        xdatai = [ x + offset*i for x in xdata ]
        ydata = ydataList[ref]
        if scale > 0:
            ydata = [ float(y)/10**scale for y in ydata ]
        for j, x in enumerate(xdatai):
            binXdataList[j].append(x)
            binYdataList[j].append( ydata[j] )
    for i in xrange( len(binXdataList) ):
        axes.plot( binXdataList[i], binYdataList[i], color="#CCCCCC", linestyle='-', linewidth=0.005 )
    
    #Draw main plots:
    for i, ref in enumerate(refs):
        xdatai = [ x + offset*i for x in xdata ]
        ydata = ydataList[ref]
        if scale > 0:
            ydata = [ float(y)/10**scale for y in ydata ]
        
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
    
    #Draw vertical lines to separate each sample:
    #for i in xrange(1, len(samples)):
    #    d = (1 - offset*len(refs))/2.0
    #    x = [i - d, i - d]
    #    y = [miny , maxy]
    #    axes.plot(x,y, color="#CCCCCC", linewidth=0.005)
    
    axes.set_xlim(xmin, xmax )
    axes.set_ylim( miny, maxy )
    libplot.editSpine( axes )

    axes.set_xticks( [ i + offset*(len(refs)/2.0) for i in range(0, len(samples))] )
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
    axes.set_ylabel( 'Percentage of mapping difference between C. Ref. and %s' % otherRefName)
    if scale > 0:
        axes.set_ylabel( 'Event counts (x%d)' %(10**scale) )
    #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )

#============draw catusRef2:
def drawRef2(rexps, exps, options, outfile, numCats):
    options.out = outfile
    fig, pdf = libplot.initImage( 8.0, 10.0, options )
    axes = fig.add_axes( [0.12, 0.14, 0.85, 0.8] )

    if len(rexps) < 1:
        return
    
    sampleNotherRefmapped = []
    ref = ''
    for sample in rexps:
        if sample == 'average':
            continue
        e = rexps[sample]
        ref = e.ref
        sampleNotherRefmapped.append( (sample, e.total) )

    otherRefName = libplot.properName( ref )
    #Set title:
    #axes.set_title("Mapability of C. Ref. in Comparison to %s" % otherRefName)
    #HACK
    axes.set_title("Mapability of C. Ref. in Comparison to GRCh37 haplotypes")
    sampleNotherRefmapped = sorted( sampleNotherRefmapped, key=lambda item: item[1], reverse=True )
    samples = [ item[0] for item in sampleNotherRefmapped]
    samples.append( 'average' )

    xdata = range( 0, len(samples) )
    colors = libplot.getColors4()
    c = -1
    #c = 0
    lines = []
    #titleDict = {'mapped':'Mapped', 'uniquelyMapped':'Uniquely Mapped', 'properlyPaired':'Properly Paired', 'uniquelyMappedAndProperlyPaired':'Uniquely Mapped And Properly Paired', 'snps':'Snp'}
    titleDict = {'mapped':'Mapped', 'properlyPaired':'Properly Paired', 'uniquelyMapped':'Uniquely Mapped', 'uniquelyMappedAndProperlyPaired':'Uniquely Mapped And Properly Paired'}
    ydataList, miny, maxy = getData2(samples, rexps, exps, titleDict.keys())
    #ydataList, miny, maxy = getData2(samples, exps, titleDict.keys())
    
    #refs = sorted( ydataList.keys() )
    offset = 0.12
    scale = -1
    if miny > 1000:
        scale = len( str(int(miny)) ) - 1

    linenames = []
    categories = ["mapped", "properlyPaired", "uniquelyMapped", "uniquelyMappedAndProperlyPaired"]
    cats = categories[:numCats]
    for i, key in enumerate( cats ):
        xdatai = [ x + offset*i for x in xdata ]
        ydata = ydataList[key]
        if scale > 0:
            ydata = [ float(y)/10**scale for y in ydata ]
        
        c += 1
        l = axes.plot( xdatai, ydata, color=colors[c], marker='.', markersize=16.0, linestyle='none')
        lines.append(l)
        linenames.append( titleDict[key] )

    if scale > 0:
        miny = float(miny)/10**scale
        maxy = float(maxy)/10**scale

    #Draw horizontal line at y = 0:
    xmin = -0.4
    xmax = len(samples) - 1 + offset*len(linenames) + offset
    axes.plot( [xmin, xmax], [0,0], color="#6B6B6B", linewidth=0.005)

    fontP = FontProperties()
    fontP.set_size('x-small')
    
    yrange = maxy - miny
    miny = miny - yrange*0.05
    maxy = maxy + yrange*0.2
    
    #Draw vertical lines to separate each sample:
    for i in xrange(1, len(samples)):
        d = (1 - offset*len(linenames))/2.0
        x = [i - d, i - d]
        y = [miny , maxy]
        axes.plot(x,y, color="#CCCCCC", linewidth=0.005)
    
    axes.set_xlim(xmin, xmax )
    axes.set_ylim( miny, maxy )
    #HACK:
    #axes.set_ylim( -2, 0 )
    libplot.editSpine( axes )

    axes.set_xticks( [ i + offset*(len(linenames)/2.0) for i in range(0, len(samples))] )
    axes.set_xticklabels( samples )
    for label in axes.xaxis.get_ticklabels():
        label.set_rotation(90)
    axes.xaxis.set_ticks_position( 'bottom' )
    axes.yaxis.set_ticks_position( 'left' )
    
    legend = pyplot.legend( lines, linenames, numpoints=1, loc='upper right', prop=fontP)
    legend._drawFrame = False

    axes.set_xlabel( 'Samples' )
    axes.set_ylabel( 'Percentage of mapping difference between C. Ref. and %s' % otherRefName) #NEED TO DO
    #axes.set_ylabel( 'Percentage of mapping difference between C. Ref. and GRCh37 haplotypes')
    if scale > 0:
        axes.set_ylabel( 'Event counts (x%d)' %(10**scale) )
    #axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    libplot.writeImage( fig, pdf, options )


def drawPlots(options, outdir, exps, rexps):
    mappedFile = os.path.join(outdir, 'mappingStats-mapped')
    drawPlot(rexps, exps, options, mappedFile, "mapped")

    umappedFile = os.path.join(outdir, "mappingStats-uniqMapped")
    drawPlot(rexps, exps, options, umappedFile, "uniquelyMapped")

    ppairedFile = os.path.join(outdir, "mappingStats-properlyPaired")
    drawPlot(rexps, exps, options, ppairedFile, 'properlyPaired')

    uppairedFile = os.path.join(outdir, "mappingStats-uniqProperlyPaired")
    drawPlot(rexps, exps, options, uppairedFile, 'uniquelyMappedAndProperlyPaired')

    snpFile = os.path.join(outdir, "mappingStats-snp")
    drawPlot(rexps, exps, options, snpFile, 'snps')

    file = os.path.join(outdir, "mappingStats-snpSingle")
    drawSamplePlot(rexps, exps, options, file, 'snps')

    for i in [1,2,3,4]:
        file = os.path.join(outdir, "mappingStats-consensusRef2-%d" %i)
        drawRef2(rexps, exps, options, file, i)

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

def tab( f, exps, rexps, samples ):
    for sample in samples:
        expList = copy.copy(exps[sample])
        expList.sort()
        expList.append( rexps[sample] )
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

            if e.ref != 'cactusRef':
                f.write("& %s & %s & %s & %s & %s & %s \\\\\n" %(ref, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
            elif e.ref == 'cactusRef' and e.weight == 2:
                f.write("& \\cellcolor{cyan!30} %s%d & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s & \\cellcolor{cyan!30} %s \\\\\n" %(ref, e.weight, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
            else:
                f.write("& %s%d & %s & %s & %s & %s & %s \\\\\n" %(ref, e.weight, libplot.prettyInt(e.mapped), libplot.prettyInt(e.uniquelyMapped), libplot.prettyInt(e.properlyPaired), libplot.prettyInt(e.uniquelyMappedAndProperlyPaired), libplot.prettyInt(e.snps)))
        f.write("\\hline\n")

def makeLatexTab( outfile, exps, rexps ):
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

    tab( f, exps, rexps, samples )
    captionStr = "" 
    label = ""
    libplot.tableCloser(f, captionStr, label)
    libplot.writeDocumentEnd( f )
    f.close()

#================= END OF MAKING LATEX TABLE ===============

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Input directory. Required argument.')
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory. Default = ./mapStats')
    #parser.add_option('-r', '--refs', dest='refs', help='Comma separated list of refs', default='hg19,apd,cox,dbb,mann,mcf,qbl,ssto,venter')
    #parser.add_option()

def checkOptions( args, options, parser ):
    if options.indir == None:
        parser.error('No input directory was given.\n')
    if not os.path.exists( options.indir ):
        parser.error('Input directory does not exist.\n')
    if options.outdir == None:
        options.outdir = os.path.join( os.getcwd(), "mapStats" )
    #options.refs = options.refs.split(',')
    #system("mkdir -p %s" % options.outdir)

def main():
    usage = "mappingPlot.py [options]"
    parser = OptionParser(usage = usage)
    initOptions( parser )
    libplot.initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( args, options, parser )
    libplot.checkOptions( options, parser )

    refIndir = os.path.join(options.indir, "otherRefs")
    refExps = getRefStats( refIndir )
    #DEBUG:
    #for r in refExps:
    #    sys.stderr.write("Ref: %s\t" %r)
    #    for s in refExps[r]:
    #        sys.stderr.write("%s-%s\t"%(s, refExps[r][s].sample ))
    #    sys.stderr.write("\n")
    
    expSets, setnames = getExpSets(options.indir)
    for i, dirs in enumerate(expSets):
        name = setnames[i]
        exps = getStats(options.indir, dirs)
        #DEBUG:
        #for s in exps:
        #    sys.stderr.write("Sample %s\t" %s)
        #    for exp in exps[s]:
        #        sys.stderr.write("%s-%d\t" %(exp.ref, exp.weight))
        #    sys.stderr.write("\n")
        
        outdir = os.path.join( options.outdir, name )
        system("mkdir -p %s" %(outdir))
        for r in refExps:
            rexps = refExps[r]
            routdir = os.path.join(outdir, r)
            system("mkdir -p %s" %(routdir))
            makeLatexTab( os.path.join(routdir, "mappingStats.tex"), exps, rexps )
            drawPlots(options, routdir, exps, rexps)

if __name__ == '__main__':
    main()
