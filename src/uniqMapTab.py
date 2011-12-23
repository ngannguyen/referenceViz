#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
import libPlotting as libplot
from sonLib.bioio import system

#Input directory has the structure:
#indir/
#   sample/
#       ref/
#           hg19/
#               lenDist.txt
#               repeat.txt
#               uniqMapStats.txt
#           cref/
#               similary structure to hg19/

def tabHeader(f):
    f.write("\\begin{table}\n") 
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|r|r|r|r|r}\n")
    #f.write("\\multicolumn{6}{c}{%s} \\\\\n" %title)
    #f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("Sample & Reads & Total Bases & \\% Repeats & SNP Rate & Enriched SNP Ratio\\\\\n")
    #f.write("Sample & Reads & Total Bases & \\%%Repeats & SNP Rate & Overal Snp Rate\\\\\n")
    f.write("\\hline\n")

def tab(f, data ):
    altColor = 1
    sortedsamples = sorted( data.keys() )
    for i,s in enumerate(sortedsamples):
        if s == 'average':
            sortedsamples.pop(i)
    sortedsamples.append('average')

    for sample in sorted( data.keys() ):
        stats = data[sample]
        if altColor == 1:
            f.write( "%s & %s & %s & %.2f \\%% & %.4f & %.4f\\\\\n" %(sample, libplot.prettyInt(stats[0]), libplot.prettyInt(stats[1]), stats[2], stats[3], stats[3]/stats[4]) )
        else:
            f.write( "\\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %.2f \\%% & \\cellcolor[gray]{0.9} %.4f & \\cellcolor[gray]{0.9} %.4f\\\\\n" %(sample, libplot.prettyInt(stats[0]), libplot.prettyInt(stats[1]), stats[2], stats[3], stats[3]/stats[4]) )
        altColor = 1 - altColor
    f.write("\\hline\n")
    return

def makeLatexTab(data, outfile):
    #data{ sample: [numreads, ubases, repeat, snprate, totalsnprate] }
    f = open(outfile, 'w')
    libplot.writeDocumentStart(f)
    #title = "Some title"
    tabHeader(f)
    #tabHeader(f, title)
    tab( f, data )
    captionStr = ""
    label = ""
    libplot.tableCloser(f, captionStr, label)
    libplot.writeDocumentEnd(f)
    f.close()

def readRepeat(file):
    if not os.path.exists(file):
        sys.stderr.write("File %s does not exist!\n" %file)
        sys.exit(1)
    f = open(file, 'r')
    items = f.readline().split()
    if len(items) < 3:
        sys.stderr.write("File %s has wrong format. Required 3 fields\n" %(file))
        sys.exit(1)
    f.close()
    return float(items[2])

def readStatsFile(file):
    if not os.path.exists(file):
        sys.stderr.write("File %s does not exists!\n" %file)
        sys.exit(1)
    f = open(file, 'r')
    numreads = int( f.readline().strip().split(': ')[1] )
    f.readline()
    #Total
    totalItems = f.readline().split('\t')
    totalbases = int(totalItems[2])
    totalSnpRate = 0
    if totalbases > 0:
        totalSnpRate = float(totalItems[1])/totalbases
    #Uniq
    items = f.readline().split('\t')
    ubases = int(items[2])
    uSnpRate = 0
    if ubases > 0:
        uSnpRate = float(items[1])/ubases
    f.close()
    return numreads, ubases, uSnpRate, totalSnpRate

def getData(indir):
    stats = {} #stats{ ref:{ hg19/cref: { sample: [numreads,ubases,repeat,usnprate,totalsnprate]} } }
    samples = os.listdir(indir)
    cref = 'cref'
    for sample in samples:
        sampledir = os.path.join(indir, sample)
        refs = os.listdir(sampledir) #apd, cox, hg19, ...
        for ref in refs:
            refpair = [ref, cref] #hg19, cref
            if ref not in stats:
                stats[ref] = {}
            
            refdir = os.path.join(sampledir, ref)
            for r in refpair:
                if r not in stats[ref]:
                    stats[ref][r] = {}
                dir = os.path.join(refdir, r)
                repeatFile = os.path.join(dir, 'repeat.txt')
                repeat = readRepeat(repeatFile)
                statsFile = os.path.join(dir, 'uniqMapStats.txt')
                numreads, ubases, usnprate, totalsnprate = readStatsFile(statsFile)
                stats[ref][r][sample] = [numreads, ubases, repeat, usnprate, totalsnprate]
    #Average:
    for ref in stats:
        for r in stats[ref]:
            avr = [0,0,0,0,0]
            for sample in stats[ref][r]:
                for i, val in enumerate( stats[ref][r][sample] ):
                    avr[i] += val
            numsample = len( stats[ref][r].keys() )
            if numsample > 0:
                for i in xrange( len(avr) ):
                    avr[i] = avr[i]/numsample
            stats[ref][r]['average'] = avr
    return stats

def main():
    usage = 'Usage %prog [options]'
    parser = OptionParser()
    parser.add_option('-i', '--indir', dest='indir', help="Required argument. Input directory")
    parser.add_option('-o', '--outdir', dest='outdir', default = '.', help="Output directory. Default is current directory")
    options, args = parser.parse_args()

    if not options.indir:
        parser.error("Please specify input directory.\n")
    if not os.path.isdir(options.indir):
        parser.error("Input directory is not a directory: %s\n" %options.indir)
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)
    stats = getData(options.indir)
    for ref in stats:
        for r in stats[ref]:
            outfile = os.path.join(options.outdir, '%s-%s-uniqMapTab.tex' %(ref, r))
            makeLatexTab(stats[ref][r], outfile)

if __name__ == '__main__':
    main()

