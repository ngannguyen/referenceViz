#!/usr/bin/env python

'''
Thu Jul 11 14:52:43 PDT 2013
nknguyen soe ucsc edu
Wrapper to get geneCompare stats for all pairs of samples.

Input: <sampleList> <genelistDirectory> <annotatedAlignmentDirectory> <geneBedDir> <halFile> <outputDirectory>
Where the formats are as following:
    genelistDirectory/
        sample1.lst (format: <geneName>\t<geneLength>)
        sample2.lst
        ...
    annotatedAlignmentDirectory/
        target1/
            query1/
                gene1.psl, gene2.psl, etc...
            query2/
                ...
            ...
        target2/
            ...
        ...
    
    geneBedDir/ (bed files contain gene annotation of each sample)
        sample1.bed
        sample2.bed
        ...
Output:
    outputDirectory/
        target1/
            query1/
            ...
        target2/
        ...
'''

import re, os, sys
from optparse import OptionParser
import numpy as np

from sonLib.bioio import system

def getHalAlignments(halfile, query, queryBedFile, target, outdir):
    f = open(queryBedFile, 'r')
    for line in f: #each gene
        line = line.strip()
        items = line.split('\t')
        assert len(items) >= 4
        gene = items[3]
        
        tempBedFile = os.path.join(outdir, "%s.bed" %gene)
        tf = open(tempBedFile, 'w')
        tf.write("%s\n" %line)
        tf.close()

        #get halLiftover
        targetPslFile = os.path.join(outdir, "%s.psl" %gene)
        cmd = 'halLiftover %s %s %s %s %s --outPSL' %(halfile, query, tempBedFile, target, targetPslFile)
        system(cmd)
        system("rm -R %s" %tempBedFile)
    f.close()

def getStats(samples, genelistdir, alndir, beddir, halfile, outdir, options):
    halAlnDir = os.path.join(outdir, "halPsls")
    system("mkdir -p %s" %halAlnDir)

    for target in samples:
        for query in samples:
            if target != query:
                #get halAlignments:
                queryBedFile = os.path.join(beddir, "%s.bed" %query)
                currHalAlnDir = os.path.join(halAlnDir, target, query)
                system("mkdir -p %s" %currHalAlnDir)
                getHalAlignments(halfile, query, queryBedFile, target, currHalAlnDir)
                
                #compare hal alignments with the annotated alignments: 
                genelistfile = os.path.join(genelistdir, query)
                currAlndir = os.path.join(alndir, target, query)
                currOutdir = os.path.join(outdir, target, query)
                system("mkdir -p %s" %currOutdir)
                cmd = "geneCompare.py %s %s %s -o %s --minCoverage %f --minAgreement %f" \
                     %(genelistfile, currAlndir, currHalAlnDir, currOutdir, options.minCoverage, options.minAgreement)
                system(cmd)
        
        #aggregateTargetStats()
    #aggregateAllStats()
                
def readList(file):
    l = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) > 0 and line[0] != '#':
            l.append(line)
    f.close()
    return l

def addOptions(parser):
    parser.add_option('--minCoverage', dest='minCoverage', type='float', default=0.9,\
                       help='Minimum portion of the query gene mapped to the target\
                       to be considered "well" mapped. Default = %default')
    parser.add_option('--minAgreement', dest='minAgreement', type='float', default=0.9,\
                       help='Minimum agreement between the Aligner and the annotated \
                       alignment. Default = %default')

def checkOptions(parser, options, args):
    if len(args) < 6:
        parser.error("Require 6 input arguments. Only found %d.\n%s\n" % (len(args), parser.get_usage()))

def main():
    usage = "Usage: %prog [options] <sampleList> <genelistDir> <annotatedAlignmentDir> <geneBedDir> <halFile> <outDir>"
    parser = OptionParser(usage=usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, options, args)
    
    samples = readList(args[0])
    genelistdir = args[1]
    alndir = args[2]
    beddir = args[3]
    halfile = args[4]
    outdir = args[5]

    getStats(samples, genelistdir, alndir, beddir, halfile, outdir, options)

if __name__ == '__main__':
    main()


