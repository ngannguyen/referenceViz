#!/usr/bin/env python

#nknguyen soe ucsc edu
#March 26 2012
#From VCF file, contruct contigs and alignment constraints that represent the variations
#
#Input: Vcf file
#Output: Fasta file containing 1 sequence per variation
#        Cigar file containing alignment constraints


import os, sys, re, time, gzip, copy
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger

######################################
#======= MAIN PIPELINE ==============#
######################################
class Setup( Target ):
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        variants = readVcfFile(self.options.vcffile)
        globalTempDir = self.getGlobalTempDir() 
        fadir = os.path.join(globalTempDir, "fastas")
        system("mkdir -p %s" %fadir)
        cigdir = os.path.join(globalTempDir, "cigars")
        system("mkdir -p %s" %cigdir)

        batchsize = 50
        for i in xrange(0, len(variants), batchsize):
            lastindex = min( [i+batchsize, len(variants)] )
            currVariants = variants[i: lastindex]

            prevVariant = None
            nextVariant = None
            if i > 0:
                prevVariant = variants[i - 1]
            if lastindex < len(variants):
                nextVariant = variants[lastindex]
            
            faOutfile = os.path.join(fadir, "%d.fa" %i)
            cigOutfile = os.path.join(cigdir, "%d.cig" %i)
            self.addChildTarget( VariantSnippet(self.options, faOutfile, cigOutfile, currVariants, prevVariant, nextVariant) )
        
        self.setFollowOnTarget( MergeFiles(self.options.outdir, fadir, cigdir) )            
        
class VariantSnippet( Target ):
    '''For each variant, extract the 'contig' that contains that variant, and build the alignment constraint of the left and the right side of the variant.
    More specifically: Go to the left 50 bases (or 200 if it lies withing the repetitive region) or until hit the previous variant. Extract this left sequence from the reference sequence.
                       Do same thing with the right.
                       variantSnippet (contig) =  leftSeq + variant + rightSeq 
    Print cigar alignment of leftSeq-leftSeq, rightSeq-rightSeq as constraint alignments.
    '''
    def __init__(self, options, faOutfile, cigOutfile, variants, prevVariant, nextVariant):
        Target.__init__(self)
        self.options = options
        self.faOutfile = faOutfile
        self.cigOutfile = cigOutfile
        self.variants  = variants
        self.prev = prevVariant
        self.next = nextVariant
    
    def run(self):
        ff = open(self.faOutfile, 'w')
        cf = open(self.cigOutfile, 'w')
        
        for i,v in enumerate(self.variants):
            if i == 0:
                prev = self.prev
            else:
                prev = self.variants[i - 1]

            if i == len(self.variants) -1 :
                next = self.next
            else:
                next = self.variants[i + 1]

            #varheader, varseq, leftcig, rightcig
            #varList = getVariantSnippet_noOverlap(v, prev, next, self.options.faheader2seq)
            varList = getVariantSnippet(v, self.options.faheader2seq)
            
            for (varheader, varseq, leftcig, rightcig) in varList:
                #write to fastafile:
                if varseq != '':
                    ff.write(">%s\n" %varheader)
                    ff.write("%s\n" %varseq)

                #write to cigar file:
                if leftcig != '':
                    cf.write("%s\n" %leftcig)
                if rightcig != '':
                    cf.write("%s\n" %rightcig)
        
        ff.close()
        cf.close()

class MergeFiles(Target):
    def __init__(self, outdir, fadir, cigdir):
        Target.__init__(self)
        self.outdir = outdir
        self.fadir = fadir
        self.cigdir = cigdir

    def run(self):
        localTempDir = self.getLocalTempDir()
        localfafile = os.path.join(localTempDir, "indelContigs.fa")
        system("cat %s/*.fa > %s" %(self.fadir, localfafile))
        fafile = os.path.join(self.outdir, "indelContigs.fa")
        filterRedundantFasta(localfafile, fafile)

        cigfile = os.path.join(self.outdir, "indelConstraints.cig")
        system("cat %s/*.cig > %s" %(self.cigdir, cigfile))

######################################
#========== STRUCTURES ==============#
######################################
class Variant():
    def __init__(self, line):
        line = line.strip()
        items = line.split("\t")
        if len(items) < 5:
            raise ValueError("Wrong vcf format. Expected 8 fields, only see %d fields. Line: %s\n" %(len(items), line))
        self.chr = items[0].strip('chr')
        self.pos = int(items[1]) -1 #convert to base 0
        self.id = items[2] #id should have the format: spc.chr
        #if self.id == '.':
        #    self.id = "%s-%d" %(self.chr, self.pos)
        #self.id = self.id.replace(".", "_")
        self.ref = items[3]
        #self.alt = items[4]
        #if self.alt == "<DEL>":
        #    self.alt = ''
        
        sep = ','
        if re.search("/", items[4]):
            sep = "/"
        alts = items[4].rstrip(sep).split(sep)
        #Remove reference allele from observed alleles:
        self.alts = []
        for a in alts:
            if a != self.ref:
                self.alts.append(a)
        
        for i, a in enumerate( self.alts ):
            if a == "<DEL>" or a == '-':
                self.alts[i] = ''

        #if len(items) > 5:
        #self.qual = items[5]
        #self.filter = items[6]
        #self.info = items[7]
        #infoitems = self.info.split(";")
        #info = {}
        #for item in infoitems:
        #    l = item.split('=')
        #    if len(l) != 2: 
        #        continue
        #        #raise ValueError("Wrong vcf format\n")
        #    info[ l[0] ] = l[1]
        #self.info = info

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        else:
            return cmp(self.pos, other.pos)

class Region():
    def __init__(self, line):
        line = line.strip()
        items = line.split("\t")
        if len(items) < 3:
            raise ValueError("Wrong repeat file format. Require: <chrom>\\t<start>\\t<end>. Got: %s\n" %line)
        self.chr = items[0].lstrip('chr')
        self.start = int(items[1])
        self.end = int(items[2])

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

######################################
#========= UTILITILES FUNCTIONS =====#
######################################
def filterRedundantFasta(infile, outfile):
    infh = open(infile, 'r')
    header2seq = {}
    header = ''
    seq = ''
    for line in infh:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        if line[0] == '>':
            if header != '' and seq != '' and header not in header2seq:
                header2seq[header] = seq
            header = line
            seq = ''
        else:
            seq += line
    if header != '' and seq != '' and header not in header2seq:
        header2seq[header] = seq
    infh.close()

    outfh = open(outfile, 'w')
    for header, seq in header2seq.iteritems():
        outfh.write("%s\n" %header)
        outfh.write("%s\n" %seq)
    outfh.close()

def getVariantSnippet(variant, seqs):
    refstart, refend, header = getSeq( variant.chr, variant.pos, seqs )
    if not header:
        raise ValueError("Could not find fasta sequence for variant %s, %d.\n" %(variant.chr, variant.pos))
    seq = seqs[header]
    
    #Left snippet:
    leftend = variant.pos
    leftstart = max( [0, refstart, leftend - 50] )
    
    if leftend <= leftstart: #overlap with previous variant - ignore the left snippet.
        raise RuntimeError("leftend <= leftstart: %d <= %d. Variant %s, %d. Ref: %s.\n" %(leftend, leftstart, variant.chr, variant.pos, header) )
    else:
        leftheader, leftseq = extractSeq(leftstart, leftend, header, seq)
        if leftseq.islower() and leftstart == leftend - 50: #the whole sequence is repetitive 
            leftstart = max( [0, refstart, leftend - 200] )
            leftheader, leftseq = extractSeq( leftstart, leftend, header, seq )
        #remove Ns:
        noNleftseq = leftseq.split('N')[-1].split('n')[-1]
        leftstart += len(leftseq) - len(noNleftseq)
        leftheader, leftseq = extractSeq( leftstart, leftend, header, seq )
        #leftseq = noNleftseq
         
    #Right snippet:
    rightstart = variant.pos + len(variant.ref)
    rightend = min( [refend, rightstart + 50] )
    
    if rightend <= rightstart: #overlap with next variant - ignore the right snippet
        raise RuntimeError("rightend <= rightstart: %d <= %d. Variant %s, %d. Ref: %s.\n" %(rightend, rightstart, variant.chr, variant.pos, header))
    else:
        rightheader, rightseq = extractSeq( rightstart, rightend, header, seq)
        if rightseq.islower() and rightend == rightstart + 50:
            rightend = min( [refend, rightstart + 200] )
            rightheader, rightseq = extractSeq( rightstart, rightend, header, seq)
        #remove Ns:
        noNrightseq = rightseq.split('N')[-1].split('n')[0]
        rightend -= len(rightseq) - len(noNrightseq)
        rightheader, rightseq = extractSeq( rightstart, rightend, header, seq)
        #rightseq = noNrightseq

    #Whole variant snippet:
    varList = [] #each item of the list include (varHeader, varSeq, leftcig, rightcig)
    for i, alt in enumerate(variant.alts):
        varSeq = "%s%s%s" %(leftseq, alt, rightseq)
        varHeader = "knownVar.%s-alt%d.%d.%d.%d.1" % (variant.id, i, len(varSeq), 0, len(varSeq))
        #cigar string:
        #Example: 
        #cigar: dbb.chr6_dbb_hap3.4610396.1832899.73.1 0 73 + cox.chr6_cox_hap2.4795371.2051315.73.1 0 73 + 7102 M 73
        leftlen = leftend - leftstart
        leftcig = ''
        if leftlen > 0:
            leftcig = "cigar: %s %d %d %s %s %d %d %s %d M %d" %(leftheader, 0, leftlen, "+", varHeader, 0, leftlen, "+", leftlen, leftlen)

        rightlen = rightend - rightstart
        rightcig = ''
        if rightlen > 0:
            rightcig = "cigar: %s %d %d %s %s %d %d %s %d M %d" %(rightheader, 0, rightlen, "+", varHeader, len(leftseq) + len(alt), len(varSeq), "+", rightlen, rightlen)
        
        varList.append( (varHeader, varSeq, leftcig, rightcig) )
        
    return varList

def getVariantSnippet_noOverlap(variant, prevVariant, nextVariant, seqs):
    #Make sure to not overlap with previous variant and next variant.
    refstart, refend, header = getSeq( variant.chr, variant.pos, seqs )
    if not header:
        raise ValueError("Could not find fasta sequence for variant %s, %d.\n" %(variant.chr, variant.pos))
    seq = seqs[header]
    
    #Left snippet:
    leftend = variant.pos
    if prevVariant:
        leftstart = max( [0, refstart, leftend - 50, prevVariant.pos + len(prevVariant.ref) ] )
    else:
        leftstart = max( [0, refstart, leftend - 50] )
    
    if leftend <= leftstart: #overlap with previous variant - ignore the left snippet.
        leftheader = ''
        leftseq = ''
        leftcig = ''
        #raise RuntimeError("leftend <= leftstart: %d <= %d. Variant %s, %d. Ref: %s. PrevVarStart: %d, PrevVarEnd: %d\n" %(leftend, leftstart, variant.chr, variant.pos, header, prevVariant.pos, prevVariant.pos + len(prevVariant.ref)))
    else:
        leftheader, leftseq = extractSeq(leftstart, leftend, header, seq)
        if leftseq.islower() and leftstart == leftend - 50: #the whole sequence is repetitive 
            if prevVariant:
                leftstart = max( [0, refstart, leftend - 200, prevVariant.pos + len(prevVariant.ref)] )
            else:
                leftstart = max( [0, refstart, leftend - 200] )
            leftheader, leftseq = extractSeq( leftstart, leftend, header, seq )
         
    #Right snippet:
    rightstart = variant.pos + len(variant.ref)
    if nextVariant:
        rightend = min( [refend, rightstart + 50, nextVariant.pos] )
    else:
        rightend = min( [refend, rightstart + 50] )
    
    if rightend <= rightstart: #overlap with next variant - ignore the right snippet
        rightheader = ''
        rightseq = ''
        rightcig = ''
        #raise RuntimeError("rightend <= rightstart: %d <= %d. Variant %s, %d. Ref: %s. NextVarStart: %d, NextVarEnd: %d\n" %(rightend, rightstart, variant.chr, variant.pos, header, nextVariant.pos, nextVariant.pos + len(nextVariant.ref)))
    else:
        rightheader, rightseq = extractSeq( rightstart, rightend, header, seq)
        if rightseq.islower() and rightend == rightstart + 50:
            if nextVariant:
                rightend = min( [refend, rightstart + 200, nextVariant.pos] )
            else:
                rightend = min( [refend, rightstart + 200] )
            rightheader, rightseq = extractSeq( rightstart, rightend, header, seq)

    varList = []
    for i, alt in variant.alts:
        #Whole variant snippet:
        #varHeader = "1k.%s" %variant.id
        varSeq = "%s%s%s" %(leftseq, alt, rightseq)
        varHeader = "knownVar.%s-alt%d.%d.%d.%d.1" % (variant.id, i, len(varSeq), 0, len(varSeq))

        #cigar string:
        #Example: 
        #cigar: dbb.chr6_dbb_hap3.4610396.1832899.73.1 0 73 + cox.chr6_cox_hap2.4795371.2051315.73.1 0 73 + 7102 M 73
        if leftseq != '':
            leftlen = leftend - leftstart
            leftcig = "cigar: %s %d %d %s %s %d %d %s %d M %d" %(leftheader, 0, leftlen, "+", varHeader, 0, leftlen, "+", leftlen, leftlen)

        if rightseq != '':
            rightlen = rightend - rightstart
            rightcig = "cigar: %s %d %d %s %s %d %d %s %d M %d" %(rightheader, 0, rightlen, "+", varHeader, len(leftseq) + len(alt), len(varSeq), "+", rightlen, rightlen)
        varList.append( (varHeader, varSeq, leftcig, rightcig) )
    return varList

def getSeq(chr, pos, seqs):
    for header in seqs.keys():
        items = header.split('.')
        if len(items) != 6: #require that the header must have the format: spc.chr.chrlen.start.fraglen.strand
            continue
        currchr = items[1].lstrip('chr')
        start = int(items[3])
        end = start + int(items[4])
        if currchr == chr and start <= pos and pos < end:
            return start, end, header
    return -1, -1, None

def extractSeq(start, end, header, seq):
    #apd.chr6_apd_hap1.4622290.0.4622290.1
    items = header.split('.')
    relStart = start
    relEnd = end
    if len(items) == 6: 
        if items[5] != '1':
            raise ValueError("Require sequences on the forward strand: %s\n" %header)
        offset = int(items[3])
        relStart -= offset
        relEnd -= offset
    if relStart < 0 or relEnd < 0 or end < start:
        raise ValueError("Negative coordinates or end < start. Start = %d, end = %d\n" %(start, end))
    
    chr = header
    if len(items) > 2:
        chr = "%s.%s" % (items[0], items[1])
    chrlen = len(seq)
    if len(items) == 6:
        chrlen = items[2]
    newheader = ".".join( [ chr, chrlen, str(start), str(end - start), '1' ] )

    return newheader, seq[relStart:relEnd]

def getSeqs(indir):
    faseqs = {} #key = header, val = seq
    files = os.listdir(indir)
    for file in files:
        filepath = os.path.join(indir, file)
        seqs = readFastaFile(filepath)
        for header, seq in seqs.iteritems():
            faseqs[header] = seq
    return faseqs

def readFastaFile(file):
    seqs = {}
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        if line[0] == '>':
            if seq != '' and header != '':
                seqs[header] = seq
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if seq != '' and header != '':
        seqs[header] = seq
    f.close()
    return seqs

def readVcfFile(file):
    variants = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        v = Variant(line)
        variants.append(v)
    f.close()
    return variants

def readRepeatFile(file):
    repeats = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        r = Region(line)
        repeats.append(r)
    f.close()
    return repeats

########################################
#======== OPTIONS & MAIN ==============#
########################################
def checkOptions(parser, options, args):
    if not options.vcffile or not os.path.exists(options.vcffile):
        parser.error("VCF file is required. None was given or given file does not exists\n")
    if not options.fastadir or not os.path.isdir(options.fastadir):
        parser.error("Directory of fasta sequences is required. None was given or the given directory is not a directory.\n")
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)
    #if options.repeat and not os.path.exists(options.repeat):
    #    parser.error("Repeat file %s does not exists.\n" %options.repeat)
    #if options.repeat:
    #    options.repeats = readRepeatFile( options.repeat )
    
    options.faheader2seq = getSeqs(options.fastadir)

def addOptions(parser):
    parser.add_option('-v', '--vcffile', dest='vcffile', help='Required argument. Input VCF file.\n')
    parser.add_option('-f', '--fastadir', dest='fastadir', help='Required argument. Fasta directory.\n')
    parser.add_option('-o', '--outdir', dest='outdir', default='.',  help='Output directory. Default = %default.\n')
    #parser.add_option('-r', '--repeat', dest='repeat', help='File that specified repetitive regions masked by repeatmasker\n')
    #parser.add_option()

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)

    addOptions(parser)

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from referenceViz.src.vcf2constraints import *
    main()


