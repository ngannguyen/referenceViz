#!/usr/bin/env python

#nknguyen soe ucsc edu
#March 20 2010 (first day of spring!)
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
        self.setFollowOnTarget( BuildAlignments(gtfdir, fadir, self.options.outdir, self.options.refname2faheader) )

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
    def __init__(self, gtfdir, fadir, outdir, refname2faheader):
        Target.__init__(self)
        self.gtfdir = gtfdir
        self.fadir = fadir
        self.outdir = outdir
        self.refname2faheader = refname2faheader

    def run(self):
        #Merge the gtf exons
        gtffiles = getfiles(self.gtfdir, ".pickle")
        exons = mergeExons( [ os.path.join(self.gtfdir, f) for f in gtffiles ] )
        #Merge the fasta sequences
        fafiles = getfiles(self.fadir, ".pickle")
        seqs = mergeSeqs( [os.path.join(self.fadir, f) for f in fafiles] )

        #for each exon, extract the target sequences, and build pairwise alignment for them
        for exonid, targets in exons.iteritems():
            self.addChildTarget( BuildAlignmentsExon(exonid, targets, seqs, self.refname2faheader, self.outdir) )

        #Merge cigar files:
        #cigarfiles = getfiles(self.outdir, ".cigar")
        outfile = os.path.join( self.outdir, "all.cigar" )
        self.setFollowOnTarget(MergeFiles(self.outdir, ".cigar",  outfile))

class BuildAlignmentsExon(Target):
    def __init__(self, id, targets, seqs, refname2seqname, outdir):
        Target.__init__(self)
        self.id = id
        self.targets = targets
        self.seqs = seqs
        self.refname2seqname = refname2seqname
        self.outdir = outdir

    def run(self):
        seqs = extractTargetSeqs( self.targets, self.seqs, self.refname2seqname ) 
        #Write fasta files:
        globalTempDir = self.getGlobalTempDir()
        for header, seq in seqs.iteritems():
            file = os.path.join(globalTempDir, header)
            f = open(file, 'w')
            f.write('>%s\n' %header)
            f.write('%s\n' %seq)
            f.close()

        #Pairwise alignment
        headers = sorted( seqs.keys() )
        cigarfiles = []
        noAlignFile = os.path.join(self.outdir, "noAlign.txt")

        for i in xrange( len(headers) -1 ):
            file1 = os.path.join( globalTempDir, headers[i] ) #tempdir/h1
            file2 = os.path.join( globalTempDir, headers[i + 1] ) #tempdir/h2
            
            outfile = os.path.join(globalTempDir, "%s-%s.cigar" %(headers[i], headers[i+1])) #tempdir/h1-h2.cigar
            cigarfiles.append(outfile)
            self.addChildTarget( PairwiseAlignment(file1, file2, outfile, noAlignFile) )
        
        #for i in xrange(len(headers) - 1):
        #    file1 = os.path.join( globalTempDir, headers[i] ) #tempdir/h1
        #    for j in xrange( i+1, len(headers) ):
        #        file2 = os.path.join( globalTempDir, headers[j] ) #tempdir/h2
        #        outfile = os.path.join(globalTempDir, "%s-%s.cigar" %(headers[i], headers[j])) #tempdir/h1-h2.cigar
        #        cigarfiles.append(outfile)
        #        self.addChildTarget( PairwiseAlignment(file1, file2, outfile) )
        
        #Merge output cigar files:
        mergedOutfile = os.path.join(self.outdir, "%s.cigar" % self.id.replace(";", "_"))
        #self.setFollowOnTarget( MergeFiles(globalTempDir, ".cigar", mergedOutfile) )
        self.setFollowOnTarget( MergeFiles0(cigarfiles, mergedOutfile) )

class PairwiseAlignment(Target):
    def __init__(self, file1, file2, outfile, noAlignFile):
        Target.__init__(self)
        self.file1 = file1
        self.file2 = file2
        self.outfile = outfile
        self.noAlignFile = noAlignFile

    def run(self):
        system("lastz %s %s --output=%s --format=cigar --coverage=90" %(self.file1, self.file2, self.outfile) )
        
        #Check to make sure that output file is not empty:
        f = open(self.outfile, 'r')
        numCigars = 0
        for line in f:
            line = line.strip()
            if re.match('cigar:', line):
                numCigars += 1
        f.close()
        if numCigars == 0:
            nf = open(self.noAlignFile, 'a')
            nf.write("%s\t%s\n" %(os.path.basename(self.file1), os.path.basename(self.file2)))
            nf.close()
            #raise NoAlignmentError("%s and %s do not have any alignment that passed coverage cutoff of 90%.\n" %(file1, file2))
        elif numCigars > 1:
            mf = open( os.path.join(self.noAlignFile.rstrip( os.path.basename(self.noAlignFile) ), "multiAlign.txt") , 'a' )
            mf.write("%s\t%s\n" %(os.path.basename(self.file1), os.path.basename(self.file2)))
            mf.close()

class MergeFiles0(Target):
    def __init__(self, infiles, outfile):
        Target.__init__(self)
        self.infiles = infiles
        self.outfile = outfile

    def run(self):
        if len(self.infiles) > 0:
            system("cat %s > %s" %(" ".join(self.infiles), self.outfile) )

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

############### UTILITIES FUNCTIONS ############
def getfiles(dir, extension):
    allfiles = os.listdir(dir)
    files = []
    for file in allfiles:
        if re.search("%s$" %extension, file):
            files.append(file)
    return files

def readGtfFile(file, name2offset):
    exons = {} #key = txname + exonNumber, val = list of Exons
    f = open(file, 'r')
    for line in f:
        exon = Exon(line)
        if exon.tname in name2offset:
            offset = name2offset[ exon.tname ]
            exon.start -= offset
            exon.end -= offset
        header = "%s;%s" %(exon.txname, exon.exon)
        if header in exons:
            exons[header].append( exon )
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
    exons = {}
    for file in files:
        currExons = pickle.load( gzip.open(file, "rb") )
        for id, exonList in currExons.iteritems():
            if id in exons:
                exons[id].extend( exonList )
            else:
                exons[id] = exonList
    return exons

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

def extractTargetSeqs( targets, seqs, refname2seqname ):
    extractedSeqs = {}
    for target in targets:
        if target.tname not in refname2seqname:
            raise ValueError("Target %s does not have a specified sequence header\n" %target.tname)
        seqname = refname2seqname[ target.tname ]
        if seqname not in seqs:
            raise ValueError("Could not find sequence %s. Seqs include: %s\n" %(seqname, ','.join(seqs.keys()) ))
        header, seq = extractSeq( target.start, target.end, seqname, seqs[seqname] )
        extractedSeqs[header] = seq
    return extractedSeqs

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
    if not os.path.isdir(options.indir):
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

