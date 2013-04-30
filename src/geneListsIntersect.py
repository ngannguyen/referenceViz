#!/usr/bin/env python

#Tue Feb 12 14:19:33 PST 2013
#nknguyen soe ucsc edu
#
#Intersect multiple gene lists
#Input: directory containing gene lists (each file is the gene list of a genome/sample)
#       minimum number of samples (n) containing a gene (optional, default = all)
#Output: 1/ plot showing number of genes versus the number of samples
#        2/ list of genes shared by >= n samples (core-genome)
#        3/ list of other genes (pan-genome)

import os, sys, re, time
from optparse import OptionParser

######## ERROR CLASSES #######
class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

####### functions #########
def readGeneList(file, sample, gene2samples):
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        if line not in gene2samples:
            gene2samples[line] = [sample]
        else:
            if sample not in gene2samples[line]:
                gene2samples[line].append(sample)
    f.close()
    return

def sharedGenesDist(gene2samples, outbase, minNumSam):
    #print to output file: <Number of shared samples>\t<Number of genes>
    outfile = "%s-stats.txt" %outbase
    intersectFile = "%s-intersection" %outbase
    xoFile = "%s-xo" %outbase

    ifh = open(intersectFile, 'w')
    xofh = open(xoFile, 'w')

    numSam2numGene = {}
    for gene, samples in gene2samples.iteritems():
        numSam = len(samples)
        if numSam not in numSam2numGene:
            numSam2numGene[numSam] = 1
        else:
            numSam2numGene[numSam] += 1
        
        if numSam >= minNumSam:
            ifh.write("%s\t%s\n" %(gene, ','.join(samples)))
        else:
            xofh.write("%s\t%s\n" %(gene, ','.join(samples)))
    ifh.close()
    xofh.close()

    f = open(outfile, 'w')
    f.write("#Total %d genes\n" %len(gene2samples))
    f.write("#Samples\tCumulative Number of Genes\tNumber of Genes\n")
    cumulativeCount = 0
    for s in sorted(numSam2numGene.keys(), reverse=True):
        cumulativeCount += numSam2numGene[s]
        #f.write("%d\t%d\t%d\n" %(s, cumulativeCount, numSam2numGene[s]))
        f.write("%d\t%d\t%d\n" %(s, numSam2numGene[s], cumulativeCount))
    f.close()

####### Options ########
def addOptions(parser):
    parser.add_option('-n', '--numSamples', dest='numSam', type='int', help='Minimum number of samples required for a gene to be reported to the intersected list of genes. Default = all input samples')

def checkOptions(parser, args, options):
    if len(args) < 2:
        raise InputError("Require <inputDirectory> and <outputBasename>\n")
    if not os.path.exists(args[0]):
        raise InputError("Input directory %s does not exist.\n" %args[1])
    if not os.path.isdir(args[0]):
        raise InputError("Input directory %s is not a directory\n" %args[1])

####### main ########
def main():
    usage = '%prog <inputDirectory> <outputBasename>'
    parser = OptionParser(usage = usage)
    addOptions(parser) 
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    gene2samples = {}
    infiles = os.listdir(args[0])
    if not options.numSam:
        options.numSam = len(infiles)
    
    for infile in infiles:
        fileItems = infile.split('.')
        if len(fileItems) == 1:
            sample = infile
        else:
            sample = '.'.join( fileItems[:-1] )
        readGeneList( os.path.join(args[0], infile), sample, gene2samples )
    
    sharedGenesDist(gene2samples, args[1], options.numSam)

if __name__ == '__main__':
    from referenceViz.src.geneListsIntersect import *
    main()



















