#!/usr/bin/env python

#nguyen soe ucsc edu
#Apr 3 2012
#This script extracts exon sequences of the high-quality samples 
#(Sanger samples: pgf, apd, cox, etc) and maps them to the lower-quality samples
#to create exon constraints for these samples.
#Input: fasta files of all samples (high and low quality)
#       gtf files of the high-quality samples
#       list of target samples needed the exon constraints.
#Output: cigar files of the constraints, one for each sample
#
#Steps:
#1/ Read input files (gtf & fasta)
#2/ Extract exon sequences
#3/ Map exon sequences to each target sample.
#4/ Aggregate the cigar files


import os, sys, re, copy, time, gzip
from optparse import OptionParser
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import system
from sonLib.bioio import getTempDirectory
from sonLib.bioio import logger

############# MAIN PIPELINE ###########
class Setup(Target):
    '''Setting up the pipeline
    '''
    def __init__(self, options):
        Target.__init__(self)
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        
        #Read gtf files:
        gtfs = getfiles( self.options.gtfdir, ".gtf")
        gtfdir = os.path.join(globalTempDir, 'gtfs')
        system("mkdir -p %s" %gtfdir)
        for gtf in gtfs:
            outfile = os.path.join(gtfdir, "%s.pickle" % gtf.rstrip("gtf").rstrip(".") )
            self.addChildTarget( ReadGtf(gtf, self.options.gtfdir, outfile, self.options.refname2offset) )

        #Read query fasta files:
        qfastas = getfiles( self.options.qdir, ".fa" )
        qfadir = os.path.join(globalTempDir, 'qfastas')
        system("mkdir -p %s" %qfadir)
        for fa in qfastas:
            faoutfile = os.path.join(qfadir, "%s.pickle" %(fa.rstrip('fa').rstrip("."))) #globalTempDir/qfastas/sample.pickle
            self.addChildTarget( ReadFasta(fa, self.options.qdir, faoutfile) )

        #After finishing reading input files, extract exon sequences and build alignments:
        self.setFollowOnTarget( BuildAlignments(gtfdir, qfadir, self.options.tdir, self.options) )

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

############### BUILD ALIGNMENTS #################
class BuildAlignments(Target): 
    def __init__(self, gtfdir, qfadir, tfadir, options):
        Target.__init__(self)
        self.gtfdir = gtfdir
        self.qfadir = qfadir
        self.tfadir = tfadir
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()

        #Merge the gtf exons:
        gtffiles = getfiles(self.gtfdir, ".pickle")
        gene2ref2exons = mergeExons( [ os.path.join(self.gtfdir, f) for f in gtffiles ] )
       
        #Merge the query fasta sequences:
        qfafiles = getfiles(self.qfadir, ".pickle")
        qseqs = mergeSeqs( [os.path.join(self.qfadir, f) for f in qfafiles] )#key=header,val=seq

        #Make temporary output directory for each target sample:
        tfafiles = getfiles(self.tfadir, ".fa")
        targets = [t.rstrip("fa").rstrip(".") for t in tfafiles]
        for t in targets:
            system("mkdir -p %s" %( os.path.join(globalTempDir, t) ) ) #currTempDir/target

        #For each exon of each gene of each query, map to all the targets
        for gene, ref2exons in gene2ref2exons.iteritems():
            for ref, exonlist in ref2exons.iteritems():
                for exon in exonlist:
                    self.addChildTarget( BuildAlignmentsExon(globalTempDir, exon, ref, gene, qseqs, self.tfadir, self.options) )

        #Merge cigar files:
        self.setFollowOnTarget( MergeFiles(targets, globalTempDir, ".cig", self.options.outdir) )

class BuildAlignmentsExon(Target):
    def __init__(self, outdir, exon, ref, gene, qseqs, tfadir, options):
        Target.__init__(self)
        self.outdir = outdir
        self.exon = exon
        self.ref = ref
        self.gene = gene
        self.qseqs = qseqs
        self.tfadir = tfadir
        self.refname2seqname = options.refname2faheader
        self.coverage = options.coverage
        self.identity = options.identity
        #self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        exon = self.exon
        exonid = "Gene:%s;Tx:%s;Exon:%s" %(exon.genename, exon.txname, exon.exon)
        #Extract sequence from the query
        if exon.tname not in self.refname2seqname:
            raise ValueError("Query %s does not have a specified sequence header\n" %exon.tname)
        seqname = self.refname2seqname[ exon.tname ]
        if seqname not in self.qseqs:
            raise ValueError("Could not find sequence %s. All query seqs include: %s\n" %(seqname, ','.join(self.qseqs.keys()) ))
        header, seq = extractSeq( exon.start, exon.end, seqname, self.qseqs[seqname])
        
        #Write the exon sequence to a fasta file:
        exonFaFile = os.path.join(globalTempDir, "exon.fa")
        f = open(exonFaFile, 'w')
        f.write(">%s\n" %header)
        f.write("%s\n" %seq)
        f.close()

        #Align sequence to each target:
        tfafiles = getfiles(self.tfadir, ".fa")
        for tfile in tfafiles:
            targetFaFile = os.path.join(self.tfadir, tfile)
            targetname = tfile.rstrip("fa").rstrip(".")
            geneOutdir = os.path.join(self.outdir, targetname, self.gene)
            if not os.path.exists(geneOutdir):
                system("mkdir -p %s" %geneOutdir)
            outfile = os.path.join( geneOutdir, "%s.cig" %(header) )
            self.addChildTarget( Alignment(exonid, exonFaFile, targetFaFile, outfile, self.coverage, self.identity) )

class Alignment(Target):
    def __init__(self, id, query, target, outfile, coverage, identity):
        Target.__init__(self)
        self.id = id
        self.query = query
        self.target = target
        self.outfile = outfile
        self.coverage = coverage
        self.identity = identity
        
    def run(self):
        #Read target sequences:
        seqs = readFastaFile(self.target)
        globalTempDir = self.getGlobalTempDir()
        localOutfile = os.path.join(globalTempDir, "out.cig")
        if os.path.exists(self.outfile):
            system("rm %s" %self.outfile)
        for header, seq in seqs.iteritems():#each sequence of the target
            targetfile = os.path.join(globalTempDir, "%s.fa" %header)
            f = open(targetfile, 'w')
            f.write(">%s\n" %header)
            f.write("%s\n" %seq)
            f.close()

            outfile = os.path.join(globalTempDir, "%s.cig" %header)
            system("lastz %s[unmask] %s[unmask] --output=%s --format=cigar --ambiguous=iupac --coverage=%d --identity=%d --hspthresh=top100%%" %(targetfile, self.query, outfile, self.coverage, self.identity))
            system("cat %s >> %s" %(outfile, localOutfile))
            system("rm %s" %targetfile)
            system("rm %s" %outfile)

        #Filter cig file:
        cigars = readCigarFile(localOutfile)
        fcigs, failedCigs = filterCigars(cigars)
        bestOfFailedQueryCigs = getBestCigs(failedCigs, True) #for the exons with multiple hits, pick the best hit (best score, then most number of matches)
        fcigs.extend(bestOfFailedQueryCigs)
        printCigars(fcigs, self.outfile)

        #Exon did not map anywhere:
        if len(cigars) == 0:
            nf = open( os.path.join(self.outfile.rstrip( os.path.basename(self.outfile) ), "..", "noAlign.txt") , 'a' )
            nf.write("%s\t%s\t%s\n" %(self.id, os.path.basename(self.query), os.path.basename(self.target)))
            nf.close()
        elif len(failedCigs) > 0: #multi-mapping
            mf = open( os.path.join(self.outfile.rstrip( os.path.basename(self.outfile) ), "..", "multiAlign.txt") , 'a' )
            mf.write("\n%s\t%s\t%s\n" %(self.id, os.path.basename(self.query), os.path.basename(self.target)))
            for c in failedCigs:
                mf.write("\t%s\n" %c.line)
            mf.close()

############### MERGE FILES #######################
class MergeFiles(Target):
    def __init__(self, targets, indir, extension, outdir):
        Target.__init__(self)
        self.targets = targets
        self.indir = indir
        self.outdir = outdir
        self.extension = extension

    def run(self):
        #targets = getfiles(self.indir, "")
        for target in self.targets:
            #if target == 'gtfs' or target == 'qfastas':
            #    continue
            outdir = os.path.join(self.outdir, target)
            system("mkdir -p %s" %outdir)
            outfile = os.path.join(outdir, "%s.cig" %target)
            indir = os.path.join(self.indir, target) #tempdir/target
            genes = getfiles(indir, "") #gene directories

            for gene in genes:
                genedir = os.path.join(indir, gene)
                if not os.path.isdir(genedir):
                    continue
                infiles = getfiles(genedir, ".cig")
                geneoutfile = os.path.join(outdir, "%s.cig" %gene)
                for infile in infiles:
                    system("cat %s >> %s" %(os.path.join(genedir, infile), geneoutfile))
                system("cat %s >> %s" %(geneoutfile, outfile))

            if os.path.exists( os.path.join(indir, "noAlign.txt") ):
                system( "cp %s %s" %( os.path.join(indir, "noAlign.txt"), os.path.join(outdir, "noAlign.txt") ) )
            if os.path.exists( os.path.join(indir, "multiAlign.txt") ):
                system( "cp %s %s" %( os.path.join(indir, "multiAlign.txt"), os.path.join(outdir, "multiAlign.txt") ) )

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
        self.score = int(items[9])
        self.cigarstr = " ".join(items[10:])
        self.cat2bases = processCigStr(self.cigarstr)
        self.line = line

    def __cmp__(self, other):
        if self.score != other.score:
            return cmp(self.score, other.score)
        else:
            match = 0
            othermatch = 0
            if 'M' in self.cat2bases:
                match = self.cat2bases['M']
            if 'M' in other.cat2bases:
                othermatch = other.cat2bases['M']
            return cmp(match, othermatch)

############## UTILITIES ##############
def processCigStr(cigstr):
    cat2bases = {}
    items = cigstr.split()
    for i in xrange(0, len(items), 2):
        category = items[i]
        bases = int(items[i + 1])
        if category not in cat2bases:
            cat2bases[category] = bases
        else:
            cat2bases[category] += bases
    return cat2bases

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

def getBestCigs(cigars, ofquery):
    name2cigs = {} #if ofquery is true, key = cig.name1, val = best cigar with name1; else, key=cig.name2 and val is best cigar with that name2
    for cig in cigars:
        if ofquery:
            name = cig.name1
        else:
            name = cig.name2
        if name not in name2cigs or name2cigs[name] < cig:
            name2cigs[name] = cig
    return name2cigs.values() 

def printCigars(cigars, outfile):
    f = open(outfile, 'w')
    for cig in cigars:
        f.write("%s\n" %cig.line)
    f.close()

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
    return seqs

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

def getfiles(dir, extension):
    allfiles = os.listdir(dir)
    files = []
    for file in allfiles:
        if re.search("%s$" %extension, file):
            files.append(file)
    return files

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

################# MAIN ################
def addOptions(parser):
    parser.add_option('-o', '--outdir', dest='outdir', default='outputs', help='Output directory. Default=%default')
    parser.add_option('-g', '--gtfdir', dest='gtfdir', help='Required argument. Directory containing input gtf files.')
    parser.add_option('-q', '--querydir', dest='qdir', help='Required argument. Directory that contains fasta files of query sequences (corresponding with samples in the gtfdir).')
    parser.add_option('-t', '--targetdir', dest='tdir', help='Required argument. Directory that contains fasta files of target sequences (to map exons to).')
    parser.add_option('-r', '--refname2faHeader', dest='refname2faheader', help='Required argument. File mapping refname in the gtf file to fasta header. Format:<gtfname> <faheader>')
    parser.add_option('-c', '--offsets', dest='offsets', help='Vega uses "HGCHR6_MHC_***" (whole human chrom 6 with the *** MHC haplotype), which are different from 6_*** (just the *** MHC haplotype), this off set file specifies the offsets between 6_*** and HGCHR6_MHC_***. Format: <HGCHR6_MHC_***> <start coordinate of the MHC region>')
    parser.add_option('--coverage', dest='coverage', default=95, type='int', help='Minimum alignment coverage. Default=%default')
    parser.add_option('--identity', dest='identity', default=90, type='int', help='Minimum alignment identity. Default=%default')

def checkOptions(parser, options, args):
    #Input directory
    if not options.gtfdir:
        raise NonExistFileError("Gtf directory is required. None was given.\n")
    if not os.path.exists(options.gtfdir):
        raise NonExistFileError("Gtf directory %s does not exist.\n" %options.gtfdir)
    if not os.path.isdir(options.gtfdir):
        raise NotDirectoryError("Gtf directory %s is not a directory.\n" %options.gtfdir)
    
    #Output directory
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)
    if not os.path.isdir(options.outdir):
        raise NotDirectoryError("Output directory %s is not a directory.\n" %options.outdir)
    
    #Fasta directory
    if not os.path.exists(options.qdir):
        raise NonExistFileError("Query fasta directory %s does not exist.\n" %options.qdir)
    if not os.path.isdir(options.qdir):
        raise NotDirectoryError("Query fasta directory %s is not a directory.\n" %options.qdir)

    if not os.path.exists(options.tdir):
        raise NonExistFileError("Target fasta directory %s does not exist.\n" %options.tdir)
    if not os.path.isdir(options.tdir):
        raise NotDirectoryError("Target fasta directory %s is not a directory.\n" %options.tdir)

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

    addOptions(parser)

    options, args = parser.parse_args()
    checkOptions(parser, options, args)

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == "__main__":
    from referenceViz.src.mapGtfExons import *
    main()
