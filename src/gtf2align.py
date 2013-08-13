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

        #Create directory for "unionGenes" (representative version of the union list of genes of all samples)
        unionFa = os.path.join(self.options.outdir, "unionGenes.fa")
        unionCig = os.path.join(self.options.outdir, "unionGenes.cig")
        system("rm -Rf %s" %unionFa)
        system("rm -Rf %s" %unionCig)

        #After finishing reading input files, construct pairwise alignments:
        self.setFollowOnTarget( BuildAlignments(gtfdir, fadir, self.options.outdir, self.options.refname2faheader, self.options.coverage, self.options.identity, self.options.introns) )

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
    def __init__(self, gtfdir, fadir, outdir, refname2faheader, coverage, identity, introns):
        Target.__init__(self)
        self.gtfdir = gtfdir
        self.fadir = fadir
        self.outdir = outdir
        self.refname2faheader = refname2faheader
        self.coverage = coverage
        self.identity = identity
        self.introns = introns

    def run(self):
        #Merge the gtf exons
        gtffiles = getfiles(self.gtfdir, ".pickle")
        gene2ref2exons = mergeExons( [ os.path.join(self.gtfdir, f) for f in gtffiles ] )
        #Merge the fasta sequences
        fafiles = getfiles(self.fadir, ".pickle")
        seqs = mergeSeqs( [os.path.join(self.fadir, f) for f in fafiles] )

        #for each exon, extract the target sequences, and build pairwise alignment for them
        for genename, ref2exons in gene2ref2exons.iteritems():
            self.addChildTarget( BuildAlignmentsGene(genename, ref2exons, seqs, self.refname2faheader, self.outdir, self.coverage, self.identity, self.introns) )

        #Merge cigar files:
        #cigarfiles = getfiles(self.outdir, ".cigar")
        refs = sorted( [gtffile.rstrip('pickle').rstrip('.') for gtffile in gtffiles] )
        genenames = gene2ref2exons.keys()
        outdir = os.path.join( self.outdir, "all" )
        system("mkdir -p %s" % outdir)
        self.setFollowOnTarget(MergeFiles(self.outdir, refs, genenames, "cig",  outdir))


class BuildAlignmentsGene(Target):
    def __init__(self, genename, ref2exons, seqs, refname2seqname, outdir, coverage, identity, introns):
        Target.__init__(self)
        self.id = genename
        self.ref2exons = ref2exons
        self.seqs = seqs
        self.refname2seqname = refname2seqname
        self.outdir = outdir
        self.coverage = coverage
        self.identity = identity
        self.introns = introns

    def run(self):
        #EXON SEQUENCES
        ref2header2seq_exon = extractTargetSeqs( self.ref2exons, self.seqs, self.refname2seqname, False ) 
        #Write fasta files:
        globalTempDir = self.getGlobalTempDir()
        for ref in ref2header2seq_exon: #each reference
            system("mkdir -p %s" %( os.path.join(globalTempDir, ref)))
            for header, seq in ref2header2seq_exon[ref].iteritems(): #write each exon to a fasta file (because lastz's [multiple] option omitted extremely short alignments and therefore miss some exon. So to have each exon in a separate file, we don't have to use the [multiple] option anymore
                file = os.path.join(globalTempDir, ref, "%s__EXON" %header)
                f = open(file, 'w')
                f.write('>%s\n' %header)
                f.write('%s\n' %seq)
                f.close()

        #INTRON SEQUENCES
        if self.introns:
            ref2header2seq_intron = extractTargetSeqs( self.ref2exons, self.seqs, self.refname2seqname, True )
            for ref in ref2header2seq_intron: #each reference
                for header, seq in ref2header2seq_intron[ref].iteritems(): 
                    file = os.path.join(globalTempDir, ref, "%s__INTRON" %header)
                    f = open(file, 'w')
                    f.write('>%s\n' %header)
                    f.write('%s\n' %seq)
                    f.close()

        #Pairwise alignment
        headers = sorted( ref2header2seq_exon.keys() ) #list of reference
        cigdir = os.path.join(globalTempDir, "cigars")
        system("mkdir -p %s" %cigdir)
        noAlignFile = os.path.join(self.outdir, "noAlign.txt")

        for i in xrange( len(headers) ):
            ref1 = headers[i]
            headers1 = ["%s__EXON" % h for h in ref2header2seq_exon[ref1].keys()] 
            if self.introns:
                headers1.extend( ["%s__INTRON" % h for h in ref2header2seq_intron[ref1].keys()] )
            for j in xrange( len(headers) ):
                if i == j:
                    continue
                ref2 = headers[j]
                headers2 = ["%s__EXON" % h for h in ref2header2seq_exon[ref2].keys()]
                if self.introns:
                    headers2.extend( ["%s__INTRON" % h for h in ref2header2seq_intron[ref2].keys()] )
                self.addChildTarget( GenePairwiseAlignment(self.id, ref1, headers1, ref2, headers2, globalTempDir, cigdir, self.coverage, self.identity, noAlignFile) )

        #Merge output cigar files:
        for i in xrange( len(headers) ):
            ref1 = headers[i]
            system( "mkdir -p %s" % os.path.join(self.outdir, ref1) )
            for j in xrange( len(headers) ):
                if i == j:
                    continue
                ref2 = headers[j]
                system( "mkdir -p %s" % os.path.join(self.outdir, ref1, ref2) )
        genename = self.id.replace(";", "_")
        ext = 'cig'
        self.setFollowOnTarget( MergeFiles0(genename, cigdir, self.outdir, ref2header2seq_exon, ext) )
        
        #Find the longest (hence representative?) version of the current gene
        geneHeader, geneSeq = getLongestTx( self.ref2exons, self.seqs, self.refname2seqname )
        geneFa = os.path.join(self.outdir, "unionGenes.fa")
        geneCig = os.path.join(self.outdir, "unionGenes.cig")
        makeGeneConstraint(self.id, geneHeader, geneSeq, geneFa, geneCig)

def makeGeneConstraint(genename, header, seq, fafile, cigfile):
    #Write sequence
    ff = open(fafile, 'a')
    ff.write(">%s\n" %genename)
    ff.write("%s\n" %seq)
    ff.close()

    #Write cigar file
    #cigar: geneName 0 geneLen + header 0 geneLen + someScore M geneLen
    l = len(seq)
    cf = open(cigfile, 'a')
    cig = "cigar: %s 0 %d + %s 0 %d + %d M %d" %(genename, l, header, l, l, l)
    cf.write("%s\n" %cig)
    cf.close()

class GenePairwiseAlignment(Target):
    def __init__(self, id, ref1, headers1, ref2, headers2, indir, outdir, coverage, identity, noAlignFile):
        Target.__init__(self)
        self.id = id
        self.ref1 = ref1
        self.headers1 = headers1
        self.ref2 = ref2
        self.headers2 = headers2
        self.indir = indir
        self.outdir = outdir
        self.coverage = coverage
        self.identity = identity
        self.noAlignFile = noAlignFile

    def run(self):
        localTempDir = self.getLocalTempDir()
        for exon1 in self.headers1:
            file1 = os.path.join(self.indir, self.ref1, exon1)
            items1 = exon1.split('__')
            assert len(items1) == 2

            cigfiles = []
            for exon2 in self.headers2:
                file2 = os.path.join(self.indir, self.ref2, exon2)
                items2 = exon2.split('__')
                assert len(items2) == 2
                if items1[1] != items2[1]: #only align exon to exon and intron to intron
                    continue

                localOutfile = os.path.join(localTempDir, "%s_%s-%s_%s-%s.cig" %(self.ref1, items1[0], self.ref2, items2[0], items1[1]))
                system("lastz %s[unmask] %s[unmask] --output=%s --format=cigar --coverage=%d --identity=%d --hspthresh=top100%%" %(file1, file2, localOutfile, self.coverage, self.identity) )
                #system("lastz %s[unmask,multiple] %s[unmask,multiple] --output=%s --format=cigar --coverage=%d --identity=%d --hspthresh=top100%%" %(self.file1, self.file2, outfile, self.coverage, self.identity) )
        
                cigars = readCigarFile(localOutfile)
                if len(cigars) > 0:
                    cigfiles.append( localOutfile )
            
            if len(cigfiles) == 1: #unique match
                system("mv %s %s" %(cigfiles[0], self.outdir))
            elif len(cigfiles) == 0: #no match
                nf = open(self.noAlignFile, 'a')
                nf.write("%s\t%s\t%s\n" %( self.id, os.path.basename(file1), ",".join(self.headers2) ))
                nf.close()
                #raise NoAlignmentError("%s and %s do not have any alignment that passed coverage cutoff of 95%. Cigar file %s\n" %(self.file1, self.file2, outfile))
            else: #multiple matches
                mf = open( os.path.join(self.noAlignFile.rstrip( os.path.basename(self.noAlignFile) ), "multiAlign.txt") , 'a' )
                mf.write("\n%s\t%s\t%s\n" %(self.id, os.path.basename(file1), ','.join([os.path.basename(cigfile) for cigfile in cigfiles]) ))
                mf.close()

def splitCigarFilesByRefPair(files, ext):
    #cox_cox.chr6_cox_hap2.4795371.1422795.276.1-pgf_hg19.chr6.171115067.29911044.276.1-[EXON,INTRON].cig
    pair2files = {}
    for file in files: #ref1_exon1-ref2_exon2-EXON.cigar or ref1_exon1-ref2_exon2-INTRON.cigar
        items = [ item.split('_')[0] for item in os.path.basename(file).rstrip(ext).rstrip('.').split('-')[:-1] ]
        refs = '-'.join(items)
        if refs not in pair2files:
            pair2files[refs] = [file]
        else:
            pair2files[refs].append(file)
    return pair2files

class MergeFiles0(Target):
    def __init__(self, genename, indir, outdir, ref2header2seq, ext):
        Target.__init__(self)
        self.genename = genename
        self.indir = indir
        #self.infiles = os.listdir(indir)
        self.outdir = outdir
        self.ref2header2seq = ref2header2seq
        self.ext = ext

    def run(self):
        genename = self.genename
        self.infiles = getfiles(self.indir, self.ext)
        if len(self.infiles) == 0:
            f = open("alignStats.txt", 'a')
            f.write("%s\tNo alignment\n" % genename)
            f.close()
            return
        
        pair2files = splitCigarFilesByRefPair(self.infiles, self.ext)
        for pair, infiles in pair2files.iteritems():
            infiles = [os.path.join(self.indir, infile) for infile in infiles]
            self.addChildTarget( PairMergeFiles(genename, pair, infiles, self.outdir, self.ref2header2seq) )

class PairMergeFiles(Target):
    def __init__(self, genename, pair, infiles, outdir, ref2header2seq):
        Target.__init__(self)
        self.genename = genename
        refs = pair.split('-')
        assert len(refs) == 2
        #if len(refs) != 2:
        #    raise ValueError("pair %s, len(refs) %d, refs :%s\n" %(pair, len(refs), ';'.join(refs) ))
        self.ref1 = refs[0]
        self.ref2 = refs[1]
        self.infiles = infiles
        self.outdir = os.path.join(outdir, self.ref1, self.ref2)
        self.ref2header2seq = ref2header2seq

    def run(self):
        outfile = os.path.join(self.outdir, "%s.cig" % (self.genename))

        if os.path.exists( outfile ):
            system("rm %s" % outfile)
        for infile in self.infiles:
            system("cat %s >> %s" %(infile, outfile))

        #Filter by exon alignment coverage (NOT the whole exon + intron coverage)
        exonfile = os.path.join(self.outdir, "%s-EXON.cig" %(self.genename))
        for infile in self.infiles:
            if re.search("EXON", infile):
                system("cat %s >> %s" %(infile, exonfile))

        cigars = readCigarFile(exonfile)
        system("rm -f %s" %exonfile)

        alignedExons = []
        alignedBases = 0
        for cig in cigars:
            if cig.name1 not in alignedExons:
                alignedExons.append(cig.name1)
                alignedBases += cig.getMatchedBases()
        allexons = []
        allBases = 0
        if self.ref1 in self.ref2header2seq:
            header2seq = self.ref2header2seq[self.ref1]
            allexons = header2seq.keys()
            for header in header2seq.keys():
                r = Ref(header)
                allBases += r.fraglen 

        pc = 0.0
        if len(allexons) > 0:
            pc = 100.0*alignedBases/allBases
        if pc < 90.0:#if only < 90% of gene mapped, filter out the alignment
            system("rm %s" %outfile)
        f = open(os.path.join(self.outdir, "%s-%s_alignStats.txt") % (self.ref1, self.ref2), 'a')
        f.write("%s\t%d\t%d\t%d\t%d\t%.3f\n" %(self.genename, len(alignedExons), len(allexons), alignedBases, allBases, pc))
        f.close()

class MergeFiles(Target):
    def __init__(self, indir, refs, genes, extension, outdir):
        Target.__init__(self)
        self.indir = indir
        self.refs = refs
        self.genes = genes
        self.ext = extension
        self.outdir = outdir

    def run(self):
        for gene in self.genes:
            outfile = os.path.join(self.outdir, "%s.%s" %(gene, self.ext))
            self.addChildTarget( GeneMergeFiles(gene, self.indir, self.refs, self.ext, outfile) )
        
        allfile = os.path.join(self.outdir, "all.cig")
        self.setFollowOnTarget( MergeAll(self.outdir, self.ext, allfile) )

def getConnectedComponent(componentList, name):
    for i, comp in enumerate(componentList):
        if name in comp:
            return i, comp
    return -1, None

class GeneMergeFiles(Target):
    def __init__(self, gene, indir, refs, ext, outfile):
        Target.__init__(self)
        self.gene = gene
        self.indir = indir
        self.refs = refs
        self.ext = ext
        self.outfile = outfile

    def run(self):
        refs = self.refs
        numRef = len(refs) - 1
        visitedRefs = []
        connectedComponents = []
        for i in xrange( len(refs) - 1 ):
            if len(visitedRefs) == numRef:
                break
            ref1 = refs[i]
            for j in xrange( i+1, len(refs) ):
                if len(visitedRefs) == numRef:
                    break
                ref2 = refs[j]
                
                currdir = os.path.join(self.indir, ref1, ref2)
                genefile = os.path.join(currdir, "%s.%s" %(self.gene, self.ext))
                if os.path.exists(genefile):
                    #if ref1 and ref2 are already connected, do not add this edge
                    c1, comp1 = getConnectedComponent(connectedComponents, ref1)
                    c2, comp2 = getConnectedComponent(connectedComponents, ref2)
                    if c1 == c2:
                        if c1 == -1 : #both are new
                            connectedComponents.append([ref1, ref2])
                        else: #ref1 and ref2 are already connect
                            continue
                    else:
                        if c1 == -1: #ref1 is new
                            connectedComponents[c2].append(ref1)
                        elif c2 == -1: #ref2 is new
                            connectedComponents[c1].append(ref2)
                        else:#combine the two components
                            newcomp = comp1 + comp2
                            connectedComponents.remove(comp1)
                            connectedComponents.remove(comp2)
                            connectedComponents.append(newcomp)

                    if ref1 not in visitedRefs:
                        visitedRefs.append(ref1)
                    if ref2 not in visitedRefs:
                        visitedRefs.append(ref2)
                    system("cat %s >> %s" %(genefile, self.outfile))

class MergeAll(Target):
    def __init__(self, indir, ext, outfile):
        Target.__init__(self)
        self.indir = indir
        self.ext = ext
        self.outfile = outfile

    def run(self):
        system( "cat %s/*%s | sort -u > %s" %(self.indir, self.ext, self.outfile) )

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
        self.end = int(items[4]) #base 0 and now exclusive (instead of base 1 and inclusive)
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

    def __cmp__(self, other):
        if self.tname != other.tname:
            return cmp(self.tname, other.tname)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

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

    def getMatchedBases(self):
        items = self.cigarstr.strip().split(" ")
        matches = 0
        for i in xrange( 0, len(items) - 1, 2 ):
            category = items[i]
            count = int(items[i + 1])
            if category == 'M':
                matches += count
        return matches

class Ref():
    def __init__(self, line):
        line = line.strip()
        items = line.split('.')
        self.longname = line
        if len(items) == 6:
            self.name = '.'.join(items[0:2])
            self.chrlen = int(items[2])
            self.start = int(items[3])
            self.fraglen = int(items[4])
            self.strand = items[5]
            if self.strand != '1':
                raise ValueError("Got sequence header of negative strand: %s. Positive strand is required. Sorry.\n" %line)
        else:
            self.name = line
            self.chrlen = -1
            self.start = 0
            self.fraglen = -1
            self.strand = '.'

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

def extractTargetSeqs( ref2exons, seqs, refname2seqname, introns ):
    #extract the exon sequences (and intron sequences if introns=True)
    ref2extractedSeqs = {}
    for ref, exons in ref2exons.iteritems(): #each exon list 
        ref2extractedSeqs[ref] = {}
        for i, exon in enumerate(exons):
            if exon.tname not in refname2seqname:
                raise ValueError("Target %s does not have a specified sequence header\n" %exon.tname)
            seqname = refname2seqname[ exon.tname ]
            if seqname not in seqs:
                raise ValueError("Could not find sequence %s. Seqs include: %s\n" %(seqname, ','.join(seqs.keys()) ))
            
            if not introns: #(i.e exons)
                header, seq = extractSeq( exon.start, exon.end, seqname, seqs[seqname] )
                ref2extractedSeqs[ref][header] = seq
            elif i < len(exons) - 1: #introns
                nextexon = exons[i + 1]
                if exon.end < nextexon.start:
                    intronheader, intronseq = extractSeq(exon.end, nextexon.start, seqname, seqs[seqname])
                    ref2extractedSeqs[ref][intronheader] = intronseq
                elif nextexon.end < exon.start:
                    intronheader, intronseq = extractSeq(nextexon.end, exon.start, seqname, seqs[seqname])
                    ref2extractedSeqs[ref][intronheader] = intronseq
    return ref2extractedSeqs

def getLongestTx( ref2exons, seqs, refname2seqname ):
    currlen = 0
    s = -1
    e = -1
    seqname = ''
    
    for ref, exons in ref2exons.iteritems(): #each exon list
        sortedExons = sorted( exons )
        start = sortedExons[0].start
        end = sortedExons[-1].end
        
        #HACK, if hg19 has a copy of the gene, return that copy, otherwise return the longest copy
        if ref == 'pgf':
            s = start
            e = end
            seqname = refname2seqname[ sortedExons[0].tname ]
            break
        #END HACK

        l = end - start
        if l > currlen:
            currlen = l
            s = start
            e = end
            seqname = refname2seqname[ sortedExons[0].tname ]
    assert s > 0 and e > 0 and e > s
    header, seq = extractSeq( s, e, seqname, seqs[seqname] )
    return header, seq

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
        raise ValueError("Header: %s. Negative coordinates or end < start.\n\
                          relStart: %d, refEnd: %d, Start = %d, end = %d\n" \
                          %(header, relStart, relEnd, start, end))
    
    chr = header
    if len(items) > 2:
        chr = "%s.%s" % (items[0], items[1])
    chrlen = len(seq)
    if len(items) == 6:
        chrlen = items[2]
    newheader = ".".join( [ chr, chrlen, str(start), str(end - start), '1' ] )

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
    parser.add_option('--includeIntrons', dest='introns', action='store_true', default=False, help='If specified, also make constraints for the introns')

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

