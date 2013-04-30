#!/usr/bin/env python

'''
Fri Apr  5 09:33:35 PDT 2013
nknguyen soe ucsc edu
Searching for k-mers of the consensus reference sequence that do not map to any input sample sequence

Input: 
1/fasta directory with the structure:
    faDir/
        sample/
            seq1.fa
            seq2.fa
    
2/reference sequence name
3/k-mer size

Output:
Bed file containing k-mers that did not map to any input sequences
'''

import os, sys, re, time, random, gzip
from optparse import OptionParser
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system

################# OBJ ###################
class Sample():
    def __init__(self, name):
        self.name = name
        self.header2seq = {}
        #self.header2rcseq = {}

    def addSeq(self, header, seq):
        if header not in self.header2seq:
            self.header2seq[header] = seq
        else:
            logger.info("Warning: attempting to add an already existed sequence (%s) to sample %s\n" %(header, self.name))

    def setSeqs(self, header2seq):
        self.header2seq = header2seq
        #self.header2rcseq = getHeader2rcseq(header2seq)

################# PIPELINE #######################
class Setup(Target):
    def __init__(self, fadir, options):
        Target.__init__(self)
        self.fadir = fadir
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        seqDir = os.path.join(globalTempDir, "sequences")
        system("mkdir -p %s" %seqDir)

        #Read fasta sequences
        samples = os.listdir(self.fadir)
        for sample in samples:
            indir = os.path.join(self.fadir, sample)
            self.addChildTarget( ReadSampleSeqs(sample, indir, seqDir) )

        self.setFollowOnTarget( MakeKmers(self.options, seqDir) )

class ReadSampleSeqs(Target):
    def __init__(self, sampleName, indir, outdir):
        Target.__init__(self)
        self.indir = indir
        self.outdir = outdir
        self.sampleName = sampleName

    def run(self):
        sample = Sample(self.sampleName)
        seqfiles = [os.path.join(self.indir, file) for file in os.listdir(self.indir)]
        header2seq = readFastas(seqfiles)
        sample.setSeqs(header2seq)
        #Pickle sample:
        pickleFile = os.path.join(self.outdir, "%s" %self.sampleName)
        pickle.dump( sample, gzip.open(pickleFile, "wb") )

class MakeKmers(Target):
    def __init__(self, options, sampledir):
        Target.__init__(self)
        self.options = options
        self.sampledir = sampledir
    
    def run(self):
        globalTempDir = self.getGlobalTempDir()
        beddir = os.path.join(globalTempDir, "beds")
        system("mkdir -p %s" %beddir)

        #Get reference sequences and create kmers:
        refname = self.options.ref
        refPickleFile = os.path.join(self.sampledir, refname)
        ref = pickle.load(gzip.open(refPickleFile, "rb"))
        k = self.options.kmer
        for header, seq in ref.header2seq.iteritems():
            numkmers = len(seq) - k + 1 #n - k + 1
            #batches of 10000 kmers at a time
            batchsize = 10000
            for i in xrange(0, numkmers, batchsize):
                subseq = seq[i: min(numkmers + k -1, i + batchsize + k -1)]
                outbedfile = os.path.join(beddir, "%s-%d.bed" %(header, i))
                self.addChildTarget( CheckKmers(header, i, subseq, k, self.sampledir, refname, outbedfile) )
        
        self.setFollowOnTarget( MergeBeds(beddir, self.options.outfile) )

class CheckKmers(Target):
    def __init__(self, chr, offset, seq, k, sampledir, refname, outfile):
        Target.__init__(self)
        self.chr = chr #chromosome/ seq name of the reference
        self.offset = offset #location of current batch in relative to "chr"
        self.seq = seq
        self.k = k
        self.sampledir = sampledir
        self.refname = refname
        self.outfile = outfile

    def run(self):
        #Load sample sequences:
        samples = []
        sampleNames = os.listdir(self.sampledir)
        for s in sampleNames:
            if s != self.refname:
                sample = pickle.load( gzip.open(os.path.join(self.sampledir, s), "rb") )
                samples.append(sample)

        #For each kmer, print to bed file if it does not exist in any input sample:
        orphanKmers = [] 
        for i in xrange(0, len(self.seq) - self.k + 1):
            kmer = self.seq[i: i+ self.k]
            orphan = True
            for sample in samples:
                for header, seq in sample.header2seq.iteritems():
                    #rcseq = sample.header2rcseq[header]
                    #if re.search(kmer, seq) or re.search(kmer, rcseq):
                    if re.search(kmer, seq) or re.search(rc(kmer), seq):
                        orphan = False
                        break
                if not orphan:
                    break
            if orphan:
                orphanKmers.append(i)

        #Print orphan kmers to output file:
        f = open(self.outfile, "w")
        #ftest = open("TEST.bed", 'a')
        #ftest.write("Chr %s, Offset %d. Number of orphanKmers: %d\n" %(self.chr, self.offset, len(orphanKmers)))
        for start in orphanKmers:
            f.write("%s\t%d\t%d\n" %(self.chr, self.offset + start, self.offset + start + self.k))
            #ftest.write("%s\t%d\t%d\n" %(self.chr, self.offset + start, self.offset + start + self.k))
        f.close()

class MergeBeds(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile

    def run(self):
        cmd = "cat %s/*.bed > %s" %(self.indir, self.outfile)
        system(cmd)

################ UTILITY FUNCTIONS ###############
def rc(nt):
    #Return a reverse complement of the input sequence.
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    rcnt = ''
    for i in xrange( len(nt) -1, -1, -1):
        if nt[i] in complement:
            rcnt += complement[nt[i]] 
        else:
            rcnt += nt[i]
    return rcnt

def getHeader2rcseq(header2seq):
    header2rcseq = {}
    for h, s in header2seq.iteritems():
        header2rcseq[h] = rc(s)
    return header2rcseq

def readFasta(file):
    header2seq = {}
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        if line[0] == '>': #done with old sequence, start a new one
            if header != '' and seq != '':
                if header in header2seq:
                    sys.stderr.write("Warning: repeatitive records of sequence %s in file %s\n" %(header, file))
                header2seq[header] = seq
            items = line.split()
            header = items[0].lstrip('>')
            seq = ''
        else:
            seq += line
   
    #last sequence
    if header != '' and seq != '':
        if header in header2seq:
            sys.stderr.write("Warning: repeatitive records of sequence %s in file %s\n" %(header, file))
        header2seq[header] = seq
    
    f.close()
    return header2seq

def readFastas(files):
    header2seq = {}
    for file in files:
        h2s = readFasta(file)
        for h, s in h2s.iteritems():
            if h in header2seq:
                sys.stderr.write("Warning: repeatitive records of sequence %s in files %s\n" %(h,  ",".join(files)))
            header2seq[h] = s
    return header2seq

################# MAIN ###############
def main():
    usage = "%prog <fasta directory> "
    parser = OptionParser(usage = usage)
    parser.add_option('-r', '--ref', dest='ref', default='reference', help='Reference sample name. Default=%default')
    parser.add_option('-k', '--kmer', dest='kmer', type='int', default=20, help='size k of kmer. Default=%default')
    parser.add_option('-o', '--outfile', dest='outfile', default='out.bed', help='Output file name. Default=%default')

    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    
    if len(args) < 1:
        parser.error("Input fasta sequence directory is required\n")

    i = Stack( Setup(args[0], options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == '__main__':
    from referenceViz.src.crefKmer import *
    main()

