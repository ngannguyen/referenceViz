#!/usr/bin/env python

'''
nknguyen soe ucsc edu
Mon Jun  3 09:55:28 PDT 2013
Split a set of constraints into multiple subsets of non-overlapping constraints
Input: 1/ Cigar file <constraints.cig> containing the constraint alignments
       2/ Fasta file <contigs.fa> of the contigs representing the variants

Output: <constraints1.cig> <contigs1.fa>
        <constraints2.cig> <contigs2.fa>
        ...

        where <constraintsX.cig> contains the constraints for contigs in <contigsX.fa>
'''

import os, sys, re
from optparse import OptionParser

class Seq():
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        
        items = header.split(".")
        assert len(items) == 6
        subitems = items[1].split("-")
        assert len(subitems) >= 2
        self.chr = subitems[0]
        self.start = int( subitems[1] )
        length = int( items[4] )
        self.end = self.start + length

def seqCompare(seq1, seq2):
    if seq1.chr != seq2.chr:
        return cmp(seq1.chr, seq2.chr)
    elif seq1.start != seq2.start:
        return cmp(seq1.start, seq2.start)
    else:
        return cmp(seq1.end, seq2.end)

class Cigar():
    def __init__(self, line):
        line = line.strip()
        items = line.split()
        if len(items) < 12:
            raise ValueError("Wrong Cigar format. At least 12 fields are expected. Only saw %d. Line: %s\n" %(len(items), line) )
        #query
        self.qname = items[1]
        self.qstart = int(items[2])
        self.qend = int(items[3])
        self.qstrand = items[4]

        if not isinstance(self.qstart, int) or not isinstance(self.qend, int):
            raise ValueError("Wrong Cigar format. Required start and end (field 3 and 4) to be integer. Got: %s, %s\n" %(items[2], items[3]))
        #target
        self.tname = items[5]
        self.tstart = int(items[6])
        self.tend = int(items[7])
        self.tstrand = items[8]
        self.score = items[9]
        self.cigarstr = " ".join(items[10:])

    def getStr(self):
        return "cigar: %s %d %d %s %s %d %d %s %s %s" %(self.qname, self.qstart, self.qend, self.qstrand, self.tname, self.tstart, self.tend, self.tstrand, self.score, self.cigarstr)

#==========================================================#
#=============== READ INPUT FILES =========================#
#==========================================================#
def readCigarFile(file):
    cigars = {} #key = tname, value = [Cigars]
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        cigar = Cigar(line)
        if cigar.tname not in cigars:
            cigars[cigar.tname] = [cigar]
        else:
            cigars[cigar.tname].append(cigar)
    f.close()
    return cigars

#=============================================#
#============ split sequences ================#
#=============================================#
def isOverlap(seq1, seq2):
    if seq1.chr == seq2.chr and seq1.start < seq2.end and seq2.start < seq1.end:
        return True
    return False

def addSeq(seqSets, seq):
    added = False
    for seqSet in seqSets:
        if not isOverlap( seqSet[-1], seq ):
            seqSet.append(seq)
            added = True
            break
    if not added:
        seqSets.append( [seq] )

def splitSeqsToNonOverlap(seqs):
    seqSets = [] #subsets of non-overlapping sequences
    for seq in seqs:
        addSeq(seqSets, seq)
    return seqSets

def readFastaFile(file):
    seqs = []
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '' and seq != '':
                seqs.append( Seq(header, seq) )
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if seq != '' and header != '':
        assert header not in seqs
        seqs.append( Seq(header, seq) )
    f.close()
    return seqs

#=============================================#
#========= Write to output files =============#
#=============================================#
def writeOutputFiles(cigars, seqSets, outname):
    seqcount = 0
    cigcount = 0
    for i, seqSet in enumerate(seqSets):
        cigfile = "%s-%d.cig" % (outname, i+1)
        fafile = "%s-%d.fa" % (outname, i+1)
        cf = open(cigfile, 'w')
        ff = open(fafile, 'w')
        for seq in seqSet:
            #write to fasta file
            ff.write(">%s\n" % seq.header)
            ff.write("%s\n" %seq.sequence)
            seqcount += 1

            #write to cigar file
            if seq.header in cigars:
                cigs = cigars[seq.header]
                for cig in cigs:
                    cf.write("%s\n" %cig.getStr())
                    cigcount += 1
        cf.close()
        ff.close()
    print seqcount
    print cigcount

def main():
    usage = ('%prog <contigFastaFile> <cigarFile> <outBaseName>')
    parser = OptionParser( usage = usage )
    options, args = parser.parse_args()
    if len(args) < 3:
        parser.error("Required 3 inputs: configFastaFile, cigarFile, outBaseName")
    fafile = args[0]
    cigfile = args[1]
    outname = args[2]

    cigars = readCigarFile(cigfile)
    seqs = readFastaFile(fafile)
    seqs = sorted( seqs, cmp=seqCompare )
    seqSets = splitSeqsToNonOverlap(seqs)
    writeOutputFiles(cigars, seqSets, outname)

if __name__ == "__main__":
    main()




