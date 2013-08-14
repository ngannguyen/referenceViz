#!/usr/bin/env python

'''
Wed Aug 14 14:00:20 PDT 2013
Make tables for the ecoli paper
Tables: 1/ operon 2/ pairwise average summary 3/ SMALL summary table of operon & geneFamilies
'''

import os, sys, re
from optparse import OptionParser
import xml.etree.ElementTree as ET
import libPlotting as libplot

#======== OPERON STATS =========
def operonTabHeader( f, title ):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|c|c|c}\n")
    f.write("\\multicolumn{4}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("Genome & Present & Conserved & %% Conserved/Present\\\\\n")
    f.write("\\hline\n")

def operonTab( f, infile ):
    ifh = open(infile, 'r')
    avrPresent = 0.0
    avrConserved = 0.0
    
    for line in ifh:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        f.write("%s \\\\\n" %( '&'.join(items) ) )
        f.write("\\hline\n")
        if items[0] == 'Average':
            avrPresent = items[1]
            avrConserved = items[2]
    ifh.close()
    return avrPresent, avrConserved

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def makeOperonTab(infile, outfile):
    ofh = open(outfile, 'w')
    libplot.writeDocumentStart( ofh )
    title = 'Operons of Escherichia coli K12 MG1655 uid57779 mapped to other Escherichia coli and Shigella genomes'
    operonTabHeader( ofh, title )
    avrPresent, avrConserved = operonTab( ofh, infile )
    captionStr = ""
    label = ""
    tableCloser(ofh, captionStr, label))
    libplot.writeDocumentEnd(ofh)
    ofh.close()
    return avrPresent, avrConserved

#======= GENE STATS =======
def geneTabHeader( f, title ):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|c|c|c|c}\n")
    f.write("\\multicolumn{5}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("Genome & Total & Annotation & PCactus & %% PCactus/Annotation\\\\\n")
    f.write("\\hline\n")

def geneTab( f, infile ):
    ifh = open(infile, 'r')
    avrTotal = 0.0
    avrPairwise = 0.0
    avrCactus = 0.0
    
    for line in ifh:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        f.write("%s \\\\\n" %( '&'.join(items) ) )
        f.write("\\hline\n")
        if items[0] == 'Average':
            avrTotal = items[1]
            avrPairwise = items[2]
            avrCactus = items[3]
    ifh.close()
    return avrTotal, avrPairwise, avrCactus

def makeGeneTab(infile, outfile):
    ofh = open(outfile, 'w')
    libplot.writeDocumentStart( ofh )
    title = 'Orthologous genes'
    geneTabHeader( ofh, title )
    avrTotal, avrPairwise, avrCactus = geneTab( ofh, infile )
    captionStr = ""
    label = ""
    tableCloser(ofh, captionStr, label))
    libplot.writeDocumentEnd(ofh)
    ofh.close()
    return avrTotal, avrPairwise, avrCactus

######### COMBINE GENE and OPERON STATS  ########
def summaryTabHeader( f, title ):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|c|c|c|c}\n")
    f.write("\\multicolumn{5}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("Category & Total & Shared & Conserved & %% Conserved/Shared\\\\\n")
    f.write("\\hline\n")

def summaryTab( f, geneTotal, geneAnno, geneCactus, operonTotal, operonPresent, operonConserved ):
    genePc =  100.0*float(geneCactus)/float(geneAnno)
    operonPc = 100.0*float(operonConserved)/float(operonPresent)

    f.write("Gene Families & %s & %s & %s & %.2f \\\\\n" %(geneTotal, geneAnno, geneCactus, genePc) )
    f.write("\\hline\n")
    f.write("Operons & %s & %s & %s & %.2f \\\\\n" %(operonTotal, operonPresent, operonConserved, operonPc) )
    f.write("\\hline\n")

def makeSummaryTab(geneTotal, geneAnno, geneCactus, operonTotal, operonPresent, operonConserved, outfile):
    ofh = open(outfile, 'w')
    libplot.writeDocumentStart( ofh )
    title = 'Multiple Alignment Assessment'
    summaryTabHeader( ofh, title )
    summaryTab( ofh, geneTotal, geneAnno, geneCactus, operonTotal, operonPresent, operonConserved )
    captionStr = ""
    label = ""
    tableCloser(ofh, captionStr, label))
    libplot.writeDocumentEnd(ofh)
    ofh.close()

def main():
    usage = "%prog <operonFile> <geneFile> <outbase> <numOperons> [options]"
    parser = OptionParser(usage=usage)
    options, args = parser.parse_args()
    
    operonFile = args[0]
    geneFile = args[1]
    outbase = args[2]
    numOperons = int(args[3])
    
    avrOperonPresent, avrOperonConserved = makeOperonTab(operonFile, "%s-operon.tex" %outbase)
    avrGeneTotal, avrGenePairwise, avrGeneCactus = makeGeneTab(geneFile, "%s-gene.tex" %outbase)
    makeSummaryTab(avrGeneTotal, avrGenePairwise, avrGeneCactus, numOperons, avrOperonPresent, avrOperonConserved, "%s.tex" %outbase)

if __name__ == '__main__':
    main()
