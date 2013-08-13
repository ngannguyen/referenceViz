#!/usr/bin/env python

'''
Wed Jul 31 10:04:36 PDT 2013
nknguyen soe ucsc edu
Extract regions of the alignment that contain at least MIN number of 
input subset of genomes and at most MAX number of other genomes.
The output is in bed format.
This is useful when searching to region specific to a subclade of genomes, 
e.g a group of pathogen VS commensal genomes.

Input: 1/ maffile 2/ List of genomes 3/ outputBasename 
        4/MIN 5/ MAX 6/target genome(s) (1 output
bed file for regions of each of these genomes)
Output: Output bed files
'''

import os, sys, re, copy
from optparse import OptionParser

class Bed():
    def __init__(self, mafline, seq2genome):
        items = mafline.split('\t')
        assert items[0] == 's'
        self.chr = items[1]
        chritems = self.chr.split('.')
        if len(chritems) > 1:
            self.genome = chritems[0]
            self.chr = '.'.join(chritems[1:])
        else:
            if self.chr in seq2genome:
                self.genome = seq2genome[self.chr]
            else:
                self.genome = self.chr
            
        start = int(items[2])
        size = int(items[3])
        end = start + size
        strand = items[4]
        chrsize = int(items[5])
        if strand == '+':
            self.start = start
            self.end = end
        else:
            self.end = chrsize - start
            self.start = chrsize - end

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

def processBlockBeds(genome2beds, blockBeds, options):
    numIngroup = 0
    numOutgroup = 0
    for bed in blockBeds:
        if bed.genome in options.genomes:
            numIngroup += 1
        else:
            numOutgroup += 1
    if numIngroup >= options.minIngroup and numOutgroup <= options.maxOutgroup:
        for bed in blockBeds:
            if bed.genome in options.targetGenomes:
                genome2beds[bed.genome].append(bed)
            #else:
            #    if re.search('reference', bed.genome):
            #        print bed.genome, bed.chr, options.targetGenomes

def extract(maffile, options):
    genome2beds = {}
    #initialize genome2beds:
    for g in options.targetGenomes:
        genome2beds[g] = []

    f = open(maffile, 'r')
    
    blockBeds = []
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        if line[0] == 'a':
            if len(blockBeds) > 0:
                processBlockBeds(genome2beds, blockBeds, options)
                blockBeds = []
        elif line[0] == 's':
            bed = Bed(line, options.seq2genome)
            blockBeds.append(bed)

    if len(blockBeds) > 0:
        processBlockBeds(genome2beds, blockBeds, options)
    f.close()
    #for g, beds in genome2beds.iteritems():
    #    print g, len(beds)
    return genome2beds

def mergeBeds(beds):
    if len(beds) <=1:
        return beds
    mergedBeds = [ copy.copy(beds[0]) ]
    for bed in beds[1:]:
        currbed = mergedBeds[-1]
        if bed.chr == currbed.chr and bed.start == currbed.end:
            currbed.end = bed.end
        else:
            mergedBeds.append( copy.copy(bed) )
    return mergedBeds

def printBeds(genome2beds, outbase, color):
    for genome, beds in genome2beds.iteritems():
        if len(beds) == 0:
            print "Genome %s does not have any beds" %genome
            continue
        beds.sort()
        beds = mergeBeds(beds)
        outfile = "%s%s.bed" %(outbase, genome)
        f = open(outfile, 'w')
        for bed in beds:
            f.write("%s\t%d\t%d\t.\t%d\t.\t%d\t%d\t%s\n" %(bed.chr, bed.start, bed.end, bed.end - bed.start, bed.start, bed.end, color))
        f.close()

def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        items.append(line.strip())
    f.close()
    return items

def readMapFile(file):
    g2s = {}
    s2g = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split('\t')
        assert len(items) == 2
        g = items[0]
        seqs = items[1].split(',')
        assert g not in g2s
        g2s[g] = seqs
        for s in seqs:
            assert s not in s2g
            s2g[s] = g
    f.close()
    return g2s, s2g

def addOptions(parser):
    parser.add_option('--minIngroup', dest='minIngroup', type='int', help='Minimum number of listed genomes an aligned block must include to be printed out. Default=all listed genomes')
    parser.add_option('--maxOutgroup', dest='maxOutgroup', type='int', default=0, help='Maximum number of genomes Not listed for an aligned block to be printed out. Default=%default')
    parser.add_option('--targetGenomeFile', dest='targetGenomeFile', help='File containing list of genomes to get beds for. Default=all listed genomes')
    parser.add_option('--genome2seqs', dest='genome2seqs', help='Format: <genome>\t<comma,sep,list,of,sequences/chroms>. Need to specify this if maf sequence name does not contain the genome info')
    parser.add_option('--color', dest='color', default='0,0,0', help='Color of the bed entries. Default=%default')

def checkOptions(parser, args, options):
    if len(args) < 3:
        parser.error("Required 3 input arguments, only %d given.\n" %(len(args)))
    if not os.path.exists(args[0]):
        parser.error("Input maffile %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        parser.error("Input genomeList does not exist.\n" %args[1])
    options.genomes = readList(args[1])
    if not options.minIngroup or options.minIngroup > len(options.genomes):
        options.minIngroup = len(options.genomes)
    if options.targetGenomeFile:
        options.targetGenomes = readList(options.targetGenomeFile)
    else:
        options.targetGenomes = options.genomes
    options.seq2genome = {}
    if options.genome2seqs:
        options.genome2seqs, options.seq2genome = readMapFile(options.genome2seqs)

def main():
    usage = '%prog [options] <maffile> <genomeList> <outputBasename>' 
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)
    
    maffile = args[0]
    outbase = args[2]
    genome2beds = extract(maffile, options)
    printBeds(genome2beds, outbase, options.color)

if __name__ == '__main__':
    main()

