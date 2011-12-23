#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
import libPlotting as libplot

def tabHeader(f):
    f.write("\\begin{table}\n") 
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{l|l|l}\n")
    f.write("\\hline\n")
    f.write(" & Events/Positions & \% \\\\\n")
    f.write("\\hline\n")

def getUnmapped(stats):
    total = stats['Total'][0]
    mapped = int( stats['Mapped'][1] )
    unmapped = total - mapped
    unmappedPc = "%.3f" % (unmapped*100.0/total)
    return str(unmapped), unmappedPc

def tab(f, stats1):
    d = {'Total':'Total', 'Mapped Outside of MHC':'Mapped outside of MHC', 'Repetitive':'Repeats','Mapped To MHC':'Mapped to MHC', 'Tandem duplications':'Tandem duplications', 'Mapped':'Mapped', 'Multimapping bases':'Multimapping'}
    #rowOrders = ['Repetitive', 'Mapped To MHC', 'Mapped Outside of MHC', 'Tandem duplications', 'Multimapping bases']
    #altColor = 1
    f.write("")
    f.write("Repeats & %s & %s \\\\\n" %( stats1['Repetitive'][0], stats1['Repetitive'][1]) )
    #f.write("\\hline\n")
    f.write("Mapped to MHC & %s & %s\\\\\n" %( stats1['Mapped To MHC'][1], stats1['Mapped To MHC'][2]))
    f.write("Mapped Outside MHC & %s & %s\\\\\n" %( stats1['Mapped Outside of MHC'][1], stats1['Mapped Outside of MHC'][2]))
    #f.write("\\hline\n")
    
    unmapped1, unmappedPc1 = getUnmapped(stats1)

    f.write("Unmapped & %s & %s \\\\\n" %(unmapped1, unmappedPc1))
    f.write("Copy Number Change & NA & NA \\\\\n")
    f.write("Tandem Duplications & %s & %s\\\\\n" %(stats1['Tandem duplications'][0], stats1['Tandem duplications'][1]) )
    f.write("Multi-mapping & %s  & %s\\\\\n" %(stats1['Multimapping bases'][0], stats1['Multimapping bases'][1]) )
    #f.write("\\hline\n")
    f.write("Total & %s & %s\\\\\n" %( stats1['Total'][1], stats1['Total'][0]) )
    f.write("\\hline\n")

def makeLatexTab(stats1, outfile):
    f = open(outfile, 'w')
    libplot.writeDocumentStart(f)
    tabHeader(f)
    tab(f, stats1)
    captionStr = ""
    label = ""
    libplot.tableCloser(f, captionStr, label)
    libplot.writeDocumentEnd(f)
    f.close()

def readFile(file):
    stats = {} #key = category, val = list of stats
    f = open(file, 'r')
    line = f.readline()
    items = line.strip().split()
    totalevents = items[1].rstrip(',')
    total = int(items[2])
    stats['Total'] = [total, totalevents]
    for line in f.readlines():
        items = line.strip().split('\t')
        if len(items) <=1:
            break
        stats[ items[0] ] = [ item.rstrip('%') for item in  items[1:]]
    return stats

def main():
    file1 = sys.argv[1]
    stats1 = readFile(file1)
    outfile = 'largeIndelTab.tex'
    makeLatexTab(stats1, outfile)

if __name__ == '__main__':
    main()
