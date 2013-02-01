#!/usr/bin/env python

#nknguyen soe ucsc edu
#March 20 2012 (first day of spring!)
#Pasre GTF files and construct pairwise alignments of target regions that align to the same gene/exon
#Input: gtf files
#       fasta files of target sequences
#
#Output: pairwise alignments in the cigar format
#

import os, sys, re, time, gzip, copy
from optparse import OptionParser
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger

############### MAIN PIPELINE #################
class Setup( Target ):
    '''Setting up the pipeline
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()

        #Read the gtf files:
        gtfs = getfiles( self.options.indir, ".gtf")
        gtfdir = os.path.join(globalTempDir, 'gtfs')
        system("mkdir -p %s" %gtfdir)
        for gtf in gtfs:
            outfile = os.path.join(gtfdir, "%s.pickle" % gtf.rstrip("gtf").rstrip(".") ) #globalTempDir/gtfs/sample.pickle
            self.addChildTarget( ReadGtf(gtf, self.options.indir, outfile, self.options.refname2offset) )

        #Read fasta files:
        fastas = getfiles( self.options.fadir, ".fa" )
        fadir = os.path.join(globalTempDir, 'fastas')
        system("mkdir -p %s" %fadir)
        for fa in fastas:
            faoutfile = os.path.join(fadir, "%s.pickle" %(fa.rstrip('fa').rstrip("."))) #globalTempDir/fastas/sample.pickle
            self.addChildTarget( ReadFasta(fa, self.options.fadir, faoutfile) )

        #After finishing reading input files, construct pairwise alignments:
        self.setFollowOnTarget( BuildAlignments(gtfdir, fadir, self.options.outdir, self.options.refname2faheader, self.options.coverage, self.options.identity) )

#========= Read input files ==============
class ReadGtf(Target):
    def __init__(self, filename, indir, outfile, name2offset):
        Target.__init__(self)
        self.infile = os.path.join(indir, filename)
        self.outfile = outfile
        self.name2offset = name2offset

    def run(self):
        exons = readGtfFile( self.infile, self.name2offset )
        pickle.dump( exons, gzip.open(self.outfile, "wb") )

class ReadFasta(Target):
    def __init__(self, filename, indir, outfile):
        Target.__init__(self)
        self.infile = os.path.join(indir, filename)
        self.outfile = outfile

    def run(self):
        seqs = readFastaFile(self.infile)
        pickle.dump( seqs, gzip.open(self.outfile, "wb") )

#=========== Build Alignments ===========
class BuildAlignments(Target):
    def __init__(self, gtfdir, fadir, outdir, refname2faheader, coverage, identity):
        Target.__init__(self)
        self.gtfdir = gtfdir
        self.fadir = fadir
        self.outdir = outdir
        self.refname2faheader = refname2faheader
        self.coverage = coverage
        self.identity = identity

    def run(self):
        #Merge the gtf exons
        gtffiles = getfiles(self.gtfdir, ".pickle")
        gene2ref2exons = mergeExons( [ os.path.join(self.gtfdir, f) for f in gtffiles ] )
        #Merge the fasta sequences
        fafiles = getfiles(self.fadir, ".pickle")
        seqs = mergeSeqs( [os.path.join(self.fadir, f) for f in fafiles] )

        #for each exon, extract the target sequences, and build pairwise alignment for them
        for genename, ref2exons in gene2ref2exons.iteritems():
            self.addChildTarget( BuildAlignmentsGene(genename, ref2exons, seqs, self.refname2faheader, self.outdir, self.coverage, self.identity) )

        #Merge cigar files:
        #cigarfiles = getfiles(self.outdir, ".cigar")
        outfile = os.path.join( self.outdir, "all.cigar" )
        self.setFollowOnTarget(MergeFiles(self.outdir, ".cigar",  outfile))


class BuildAlignmentsGene(Target):
    def __init__(self, genename, ref2exons, seqs, refname2seqname, outdir, coverage, identity):
        Target.__init__(self)
        self.id = genename
        self.ref2exons = ref2exons
        self.seqs = seqs
        self.refname2seqname = refname2seqname
        self.outdir = outdir
        self.coverage = coverage
        self.identity = identity

    def run(self):
        ref2header2seq = extractTargetSeqs( self.ref2exons, self.seqs, self.refname2seqname ) 
        #Write fasta files:
        globalTempDir = self.getGlobalTempDir()
        for ref in ref2header2seq: #each reference
            system("mkdir -p %s" %( os.path.join(globalTempDir, ref)))
            for header, seq in ref2header2seq[ref].iteritems(): #write each exon to a fasta file (because lastz's [multiple] option omitted extremely short alignments and therefore miss some exon. So to have each exon in a separate file, we don't have to use the [multiple] option anymore
                file = os.path.join(globalTempDir, ref, header)
                f = open(file, 'w')
                f.write('>%s\n' %header)
                f.write('%s\n' %seq)
                f.close()

        #Pairwise alignment
        headers = sorted( ref2header2seq.keys() ) #list of reference
        cigarfiles = []
        noAlignFile = os.path.join(self.outdir, "noAlign.txt")

        for i in xrange( len(headers) -1 ):
            ref1 = headers[i]
            ref2 = headers[i + 1]
            for exon1 in ref2header2seq[ref1].keys():
                file1 = os.path.join(globalTempDir, ref1, exon1)
                for exon2 in ref2header2seq[ref2].keys():
                    file2 = os.path.join(globalTempDir, ref2, exon2)
                    outfile = os.path.join(globalTempDir, "%s_%s-%s_%s.cigar" %(ref1, exon1, ref2, exon2)) #tempdir/h1-h2.cigar
                    cigarfiles.append(outfile)
                    self.addChildTarget( PairwiseAlignment(self.id, file1, file2, outfile, noAlignFile, self.coverage, self.identity) )
        
        #Merge output cigar files:
        mergedOutfile = os.path.join(self.outdir, "%s.cigar" % self.id.replace(";", "_"))
        #self.setFollowOnTarget( MergeFiles(globalTempDir, ".cigar", mergedOutfile) )
        self.setFollowOnTarget( MergeFiles0(cigarfiles, mergedOutfile, ref2header2seq) )

#class BuildAlignmentsGene(Target):
#    def __init__(self, genename, ref2exons, seqs, refname2seqname, outdir, coverage, identity):
#        Target.__init__(self)
#        self.id = genename
#        self.ref2exons = ref2exons
#        self.seqs = seqs
#        self.refname2seqname = refname2seqname
#        self.outdir = outdir
#        self.coverage = coverage
#        self.identity = identity
#
#    def run(self):
#        ref2header2seq = extractTargetSeqs( self.ref2exons, self.seqs, self.refname2seqname ) 
#        #Write fasta files:
#        globalTempDir = self.getGlobalTempDir()
#        for ref in ref2header2seq: #each reference
#            file = os.path.join(globalTempDir, ref)
#            f = open(file, 'w')
#            for header, seq in ref2header2seq[ref].iteritems(): #each exon
#                f.write('>%s\n' %header)
#                f.write('%s\n' %seq)
#            f.close()
#
#        #Pairwise alignment
#        headers = sorted( ref2header2seq.keys() ) #list of reference
#        cigarfiles = []
#        noAlignFile = os.path.join(self.outdir, "noAlign.txt")
#
#        for i in xrange( len(headers) -1 ):
#            file1 = os.path.join( globalTempDir, headers[i] ) #tempdir/h1
#            file2 = os.path.join( globalTempDir, headers[i + 1] ) #tempdir/h2
#            
#            outfile = os.path.join(globalTempDir, "%s-%s.cigar" %(headers[i], headers[i+1])) #tempdir/h1-h2.cigar
#            cigarfiles.append(outfile)
#            self.addChildTarget( PairwiseAlignment(self.id, file1, file2, outfile, noAlignFile, self.coverage, self.identity) )
#        
#        #for i in xrange(len(headers) - 1):
#        #    file1 = os.path.join( globalTempDir, headers[i] ) #tempdir/h1
#        #    for j in xrange( i+1, len(headers) ):
#        #        file2 = os.path.join( globalTempDir, headers[j] ) #tempdir/h2
#        #        outfile = os.path.join(globalTempDir, "%s-%s.cigar" %(headers[i], headers[j])) #tempdir/h1-h2.cigar
#        #        cigarfiles.append(outfile)
#        #        self.addChildTarget( PairwiseAlignment(file1, file2, outfile) )
#        
#        #Merge output cigar files:
#        mergedOutfile = os.path.join(self.outdir, "%s.cigar" % self.id.replace(";", "_"))
#        #self.setFollowOnTarget( MergeFiles(globalTempDir, ".cigar", mergedOutfile) )
#        self.setFollowOnTarget( MergeFiles0(cigarfiles, mergedOutfile, ref2header2seq) )

class PairwiseAlignment(Target):
    def __init__(self, id, file1, file2, outfile, noAlignFile, coverage, identity):
        Target.__init__(self)
        self.id = id
        self.file1 = file1
        self.file2 = file2
        self.outfile = outfile
        self.noAlignFile = noAlignFile
        self.coverage = coverage
        self.identity = identity

    def run(self):
        localTempDir = self.getLocalTempDir()
        outfile = os.path.join(localTempDir, 'out.cig')
        system("lastz %s[unmask] %s[unmask] --output=%s --format=cigar --coverage=%d --identity=%d --hspthresh=top100%%" %(self.file1, self.file2, outfile, self.coverage, self.identity) )
        #system("lastz %s[unmask,multiple] %s[unmask,multiple] --output=%s --format=cigar --coverage=%d --identity=%d --hspthresh=top100%%" %(self.file1, self.file2, outfile, self.coverage, self.identity) )
        #system("lastz %s[multiple] %s[multiple] --output=%s --format=cigar --coverage=95" %(self.file1, self.file2, outfile) )
        
        cigars = readCigarFile(outfile)
        fcigs,failedCigs = filterCigars(cigars)
        printCigars(fcigs, self.outfile)

        #Check to make sure that output file is not empty:
        if len(cigars) == 0:
            nf = open(self.noAlignFile, 'a')
            nf.write("%s\t%s\t%s\n" %(self.id, os.path.basename(self.file1), os.path.basename(self.file2)))
            nf.close()
            #raise NoAlignmentError("%s and %s do not have any alignment that passed coverage cutoff of 90%.\n" %(file1, file2))
        elif len(cigars) > len(fcigs):
            mf = open( os.path.join(self.noAlignFile.rstrip( os.path.basename(self.noAlignFile) ), "multiAlign.txt") , 'a' )
            mf.write("\n%s\t%s\t%s\n" %(self.id, os.path.basename(self.file1), os.path.basename(self.file2)))
            for c in failedCigs:
                mf.write("\t%s\n" %c.line)
            mf.close()

class MergeFiles0(Target):
    def __init__(self, infiles, outfile, ref2header2seq):
        Target.__init__(self)
        self.infiles = infiles
        self.outfile = outfile
        self.ref2header2seq = ref2header2seq

    def run(self):
        genename = os.path.basename(self.outfile).rstrip('cigar').rstrip('.')
        
        if len(self.infiles) == 0:
            f = open("alignStats.txt", 'a')
            f.write("%s\tNo alignment\n" % genename)
            f.close()
            return
        
        #system("cat %s > %s" %(" ".join(self.infiles), self.outfile) )
        if os.path.exists( self.outfile ):
            system("rm %s" %self.outfile)
        for infile in self.infiles:
            system("cat %s >> %s" %(infile, self.outfile))

        cigars = readCigarFile(self.outfile)
        alignedExons = []
        for cig in cigars:
            if cig.name1 not in alignedExons:
                alignedExons.append(cig.name1)
            if cig.name2 not in alignedExons:
                alignedExons.append(cig.name2)
        allexons = []
        for ref, header2seq in self.ref2header2seq.iteritems():
            allexons.extend( header2seq.keys() )
        f = open("alignStats.txt", 'a')
        pc = 0.0
        if len(allexons) > 0:
            pc = 100.0*len(alignedExons)/len(allexons)
        f.write("%s\t%d\t%d\t%.3f\n" %(genename, len(alignedExons), len(allexons), pc))
        f.close()


class MergeFiles(Target):
    def __init__(self, indir, extension, outfile):
        Target.__init__(self)
        self.indir = indir
        self.ext = extension
        #infiles = getfiles(indir, extension)
        #self.infiles = [ os.path.join(indir, f) for f in infiles ]
        self.outfile = outfile

    def run(self):
        system("cat %s/*%s | sort -u > %s" %(self.indir, self.ext, self.outfile) )

############### STRUCTURES/OBJECTS ###############
class Exon():
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) != 9:
            raise GtfFormatError("Wrong Gtf format. Expected 9 fields, only see %d. Line: \n%s\n" %(len(items), line) )
        self.tname = items[0]
        #self.source = items[1]
        self.feature = items[2]
        self.start = int(items[3]) -1 #convert to base 0
        self.end = int(items[4]) #base 0 and now eexclusive (instead of base 1 and inclusive)
        #self.score = items[5]
        self.strand = items[6]
        #self.frame = int(items[7])
        genestr = items[8]
        gitems = genestr.rstrip(";").split('; ')
        iddict = {}
        for gi in gitems:
            gilist = gi.split()
            if len(gilist) != 2:
                raise GtfFormatError("Wrong Gtf format\n")
            iddict[ gilist[0] ] = gilist[1].lstrip('"').rstrip('"')
        self.txname = iddict['transcript_name']
        self.genename = iddict['gene_name']
        self.exon = iddict['exon_number']

class Cigar():
    def __init__(self, line):
        line = line.strip()
        items = line.split()
        if len(items) < 12:
            raise ValueError("Wrong Cigar format. At least 12 fields are expected. Only saw %d. Line: %s\n" %(len(items), line))
        #query
        self.name1 = items[1]
        self.start1 = int(items[2])
        self.end1 = int(items[3])
        self.strand1 = items[4]

        #target
        self.name2 = items[5]
        self.start2 = int(items[6])
        self.end2 = int(items[7])
        self.strand2 = items[8]
        self.score = items[9]
        self.cigarstr = " ".join(items[10:])
        self.line = line

############### UTILITIES FUNCTIONS ############
def readCigarFile(file):
    f = open(file, 'r')
    cigars = []
    for line in f:
        cigars.append( Cigar(line) )
    f.close()
    return cigars 

def filterCigars(cigars):
    #Filter out target (field 1) that mapped to multiple places
    filteredCigars = []
    failedCigars = []
    name1hits = {}
    name2hits = {}
    for cig in cigars:
        if cig.name1 not in name1hits:
            name1hits[cig.name1] = 1
        else:
            name1hits[cig.name1] += 1

        if cig.name2 not in name2hits:
            name2hits[cig.name2] = 1
        else:
            name2hits[cig.name2] += 1

    for cig in cigars:
        if name1hits[cig.name1] == 1 and name2hits[cig.name2] == 1:
            filteredCigars.append(cig)
        else:
            failedCigars.append(cig)
    return filteredCigars, failedCigars

def printCigars(cigars, outfile):
    f = open(outfile, 'w')
    for cig in cigars:
        f.write("%s\n" %cig.line)
    f.close()

def getfiles(dir, extension):
    allfiles = os.listdir(dir)
    files = []
    for file in allfiles:
        if re.search("%s$" %extension, file):
            files.append(file)
    return files

def addExonToList(exons, exon):
    #Make sure that the exon does not overlap with any exon in the list. If two exons are overlapped, pick the longer one
    overlapIndices = []
    for i, e in enumerate(exons):
        if e.start < exon.end and exon.start < e.end: #overlap
            overlapIndices.append(i)
    
    if len(overlapIndices) == 0:#did not have any overlapped exon
        exons.append(exon)
    elif len(overlapIndices) == 1: #overlapped with 1 exon:
        i = overlapIndices[0]
        e = exons[i]
        if exon.end - exon.start > e.end - e.start: #if exon is longer than e, replace e with exon:
            exons[i] = exon
    else: #overlapped with multiple exon, don't add it
        return

def readGtfFile(file, name2offset):
    exons = {} #key = genename, val = list of Exons
    f = open(file, 'r')
    for line in f:
        exon = Exon(line)
        if exon.tname in name2offset:
            offset = name2offset[ exon.tname ]
            exon.start -= offset
            exon.end -= offset
        #header = "%s;%s" %(exon.txname, exon.exon)
        header = exon.genename
        if header in exons:
            addExonToList( exons[header], exon )
        else:
            exons[header] = [exon]
    f.close()
    return exons

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
    #raise ValueError("Infile : %s; Seqs: %s. Len: %d\n" %(file, ','.join(seqs.keys()), len(seqs[ seqs.keys()[0] ]) ) )
    return seqs

def mergeExons(files):
    gene2ref2exons = {} #key = genename, val = {refname: list of exons}
    for file in files:#each ref
        ref = os.path.basename(file).rstrip("pickle").rstrip(".")
        currExons = pickle.load( gzip.open(file, "rb") )
        for id, exonList in currExons.iteritems(): #each gene
            if id in gene2ref2exons:
                gene2ref2exons[id][ref] = exonList
                #exons[id].append( exonList )
            else:
                gene2ref2exons[id] = { ref: exonList }
                #exons[id] = [exonList]
    return gene2ref2exons

def mergeSeqs(files):
    seqs = {}
    for file in files:
        currSeqs = pickle.load(gzip.open(file, "rb"))
        for header, seq in currSeqs.iteritems():
            seqs[header] = seq
    #raise ValueError("Input files %d, %s\n" %(len(files), ','.join(files) )) 
    #raise ValueError("There are %d fasta sequences: %s\n" %(len(seqs), ','.join( [ "%s___%d" %(h, len(s)) for h, s in  seqs.iteritems()] ) ))
    return seqs

def tname2faHeader(file): 
    tname2fa = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        if len(items) != 2:
            raise RefNameToFastaHeaderFormatError("Expected 2 fields, only got %d. Line: %s\n" %(len(items), line))
        tname2fa[items[0]] = items[1]
    f.close()
    return tname2fa

def extractTargetSeqs( ref2exons, seqs, refname2seqname ):
    ref2extractedSeqs = {}
    for ref, exons in ref2exons.iteritems(): #each exon list 
        ref2extractedSeqs[ref] = {}
        for exon in exons:
            if exon.tname not in refname2seqname:
                raise ValueError("Target %s does not have a specified sequence header\n" %exon.tname)
            seqname = refname2seqname[ exon.tname ]
            if seqname not in seqs:
                raise ValueError("Could not find sequence %s. Seqs include: %s\n" %(seqname, ','.join(seqs.keys()) ))
            header, seq = extractSeq( exon.start, exon.end, seqname, seqs[seqname] )
            ref2extractedSeqs[ref][header] = seq
    return ref2extractedSeqs

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
    
    #if newheader == 'hg19.chr6.171115067.31992900.159.1' or seq[relStart:relEnd] == '':
    #    raise ValueError("Start: %d, End: %d, relStart: %d, relEnd: %d, seqname: %s, lenSeq: %d, extracted seq: %s; seq10: %s\n" %(start, end, relStart, relEnd, newheader, len(seq), seq[relStart:relEnd], seq[:10] ))

    return newheader, seq[relStart:relEnd]

def readMhcOffsets(file):
    name2offset = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        name2offset[items[0]] = int(items[1])
    f.close()
    return name2offset

################ ERROR CLASSES ################
class NoAlignmentError(ValueError):
    pass

class RefNameToFastaHeaderFormatError(ValueError):
    pass

class NonExistFileError(Exception):
    pass

class NotDirectoryError(Exception):
    pass

class GtfFormatError(ValueError):
    pass

class MissingInputError(ValueError):
    pass

###################### MAIN #####################
def addOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', default='.', help='Input directory. Default=%default')
    parser.add_option('-o', '--outdir', dest='outdir', default='outputs', help='Output directory. Default=%default')
    parser.add_option('-f', '--fastaDir', dest='fadir', help='Required argument. Directory that contains fasta files of target sequences. Default=%default')
    parser.add_option('-r', '--refname2faHeader', dest='refname2faheader', help='Required argument. File mapping refname in the gtf file to fasta header. Format:<gtfname> <faheader>')
    parser.add_option('-c', '--offsets', dest='offsets', help='Vega uses "HGCHR6_MHC_***" (whole human chrom 6 with the *** MHC haplotype), which are different from 6_*** (just the *** MHC haplotype), this off set file specifies the offsets between 6_*** and HGCHR6_MHC_***. Format: <HGCHR6_MHC_***> <start coordinate of the MHC region>')
    parser.add_option('--coverage', dest='coverage', default=95, type='int', help='Minimum alignment coverage. Default=%default')
    parser.add_option('--identity', dest='identity', default=90, type='int', help='Minimum alignment identity. Default=%default')

def checkOptions(parser, options, args):
    #Input directory
    if not os.path.exists(options.indir):
        raise NonExistFileError("Input directory %s does not exist.\n" %options.indir)
    if not os.path.isdir(options.indir):
        raise NotDirectoryError("Input directory %s is not a directory.\n" %options.indir)
    
    #Output directory
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)
    if not os.path.isdir(options.outdir):
        raise NotDirectoryError("Output directory %s is not a directory.\n" %options.outdir)
    
    #Fasta directory
    if not os.path.exists(options.fadir):
        raise NonExistFileError("Fasta directory %s does not exist.\n" %options.fadir)
    if not os.path.isdir(options.fadir):
        raise NotDirectoryError("Fasta directory %s is not a directory.\n" %options.fadir)

    #refname2faheader
    if not options.refname2faheader:
        raise MissingInputError("Please specify the refname2faHeader file.\n")
    else:
        options.refname2faheader = tname2faHeader(options.refname2faheader)

    options.refname2offset = {}
    if options.offsets:
        options.refname2offset = readMhcOffsets(options.offsets)

def main():
    parser = OptionParser()
    Stack.addJobTreeOptions(parser)

    addOptions( parser )

    options, args = parser.parse_args()
    checkOptions( parser, options, args )

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from referenceViz.src.gtf2align import *
    main()

