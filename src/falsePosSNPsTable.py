#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
import xml.etree.ElementTree as ET
import libPlotting as libplot

class Snp():
    def __init__(self, line):
        items = line.strip().split('\t')
        self.chr = items[0]
        self.start = items[1]
        self.end = items[2]
        self.allele = items[3]
        self.refallele = items[4]
        self.annotation = items[5]
        self.items = [self.chr, self.start, self.allele, self.refallele, self.annotation]

def tabHeader( f, title ):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    f.write("\\begin{tabular}{c|c|c|c|c}\n")
    f.write("\\multicolumn{5}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("Chrom & Start & C.Ref. Allele & GRCh37 Allele & Annotation\\\\\n")
    f.write("\\hline\n")

def tab( f, snps ):
    for snp in snps:
        f.write("%s \\\\\n" %( '&'.join(snp.items) ))
    f.write("\\hline\n")

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def makeLatexTab( snps ):
    libplot.writeDocumentStart( sys.stdout )
    title = ''
    tabHeader( sys.stdout, title )
    tab( sys.stdout, snps )
    captionStr = ""
    label = ""
    tableCloser(sys.stdout, captionStr, label)
    libplot.writeDocumentEnd(sys.stdout)

def readfile():
    snps = []
    for line in sys.stdin:
        snps.append( Snp(line) )
    return snps

def main():
    #Infile: stdin, Outfile: sys.stdout
    snps = readfile()
    makeLatexTab(snps)

if __name__ == '__main__':
    main()
