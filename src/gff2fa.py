#!/usr/bin/env python

#Tue Feb 26 23:58:13 PST 2013
#nknguyen soe ucsc edu
#Extract sequences of CDSs specified in a gff file from genome sequences
#Input: gff file, genome fasta file
#Output: fasta file

import os, sys, re, time
from optparse import OptionParser

######### GFF #######
class Gff():
    '''Gff record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 9:
            raise GffFormatError("gff format requires 9 fields. Only found %d in:\n%s\n" %(len(items), line))
        self.desc = line
        self.chr = items[0] 
        self.source = items[1]
        self.feature = items[2] #E.g: CDS, gene, exon, start_codon, ...
        self.start = int(items[3]) #base 1
        self.end = int(items[4]) #inclusive
        self.score = items[5]
        self.strand = items[6]
        self.frame = items[7]
        self.group = items[8].rstrip(';')
        self.attr = {} #key = attribute name, value = attribute value
        attrItems = self.group.split(";")
        for item in attrItems:
            kv = item.split("=")
            if len(kv) != 2:
                raise GffFormatError("Field 9, wrong attribute format. Must have the format <Key=value>; read in '%s'\n" %item)
            self.attr[kv[0]] = kv[1]

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        else:
            return cmp(self.start, other.start)

class Transcript(Gff):
    '''Represent a transcript
    '''
    def __init__(self, line):
        Gff.__init__(self, line)
        self.exons = [] #list of Gff records, each represent an exon

class Gene(Gff):
    '''Represent a gene (similar to Transcript)
    '''
    def __init__(self, line):
        Gff.__init__(self, line)
        self.transcripts = [] #list of Transcripts

######## ERROR CLASSES #######
class GffFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

####### functions #########
def extractSeq(gffs, header2seq, qNameAttr):
    if len(gffs) == 0:
        raise ValueError("convertGff2psl: no gff to convert\n")
    
    if qNameAttr not in gffs[0].attr:
        raise ValueError("qNameAttr %s is not in field 9 of gff %s\n" % (qNameAttr, gffs[0].desc))
    header = gffs[0].attr[qNameAttr]
    seq = ''
    
    for gff in gffs:#each exon
        genomeseq = header2seq[ gff.chr ]
        exonseq = genomeseq[gff.start - 1: gff.end] #gff coordinate is base 1 and inclusive end
        seq += exonseq
    return header, seq

def readGff(file):
    genes = []
    f = open(file, 'r')
    for line in f:
        #Skip blank lines or comment lines:
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.strip().split('\t')
        if len(items) < 9:
	    raise GffFormatError('Wrong format:\n%s.\n Gene file must have gff \
            format, which requires 9 fields.' %line)
    
        #Gff must be sorted
        if items[2] == 'gene':
            genes.append( Gene(line) )
        elif items[2] == 'CDS':
            if len(genes[-1].transcripts) == 0:
                genes[-1].transcripts.append( Transcript(line) )
                genes[-1].transcripts[-1].exons.append( Gff(line) )
            else:
                newgff = Gff(line)
                exons = genes[-1].transcripts[-1].exons
                if newgff.strand != exons[-1].strand:
                    raise ValueError("Different strand within one transcript for gene with gff record:\n%s\n" %(newgff.desc))
                if newgff.strand == '+':
                    genes[-1].transcripts[-1].exons.append( newgff )
                else:
                    newexons = [newgff]
                    newexons.extend(exons)
                    genes[-1].transcripts[-1].exons = newexons

        #elif re.search('RNA', items[2]) or re.search('transcript', items[2]):
        #    genes[-1].transcripts.append( Transcript(line) )
        #elif items[2] == 'exon':
        #    genes[-1].transcripts[-1].exons.append( Gff(line) )

    #Return only genes with CDS(s):
    genesWithCDSs = []
    for g in genes:
        if len(g.transcripts) > 0:
            genesWithCDSs.append(g)

    f.close()
    return genesWithCDSs

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
            header = items[0].lstrip('>').split('|')[3]
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

def getGeneSeqs(genes, header2seq, outfile, qNameAttr):
    ofh = open(outfile, 'w')
    for gene in genes:
        gffs = []
        if len(gene.transcripts) == 0:
            header, seq = extractSeq([gene], header2seq, qNameAttr)
            ofh.write(">%s\n%s\n" % (header, seq) )
        else:
            for tx in gene.transcripts:
                if len(tx.exons) == 0:
                    gffs = [tx]
                else:
                    gffs = tx.exons

                if qNameAttr not in gffs[0].attr:
                    if qNameAttr not in gene.attr:
                        raise ValueError("qNameAttr %s is not in field 9 of gff %s\n" % (qNameAttr, gene.desc))
                    gffs[0].attr[ qNameAttr ] = gene.attr[ qNameAttr ]

                header, seq = extractSeq(gffs, header2seq, qNameAttr)
                ofh.write(">%s\n%s\n" % (header, seq) )
    ofh.close()

####### main ########
def main():
    usage = '%prog <gff> <genomeSeqFasta> <outputFasta>'
    parser = OptionParser(usage = usage)
    parser.add_option('--queryNameAttribute', dest='qNameAttr', default='ID', help='The attribute to extract (use as) query (gene) name from field 9 of the gff. Default = %default')
    #parser.add_option('--convertSeqName', dest='nameFile', help='Optional. File that mapped sequence (chromosome) name in the gff to a different name in the output file. Format: <seqNameInGff> <newName>')

    options, args = parser.parse_args()
    if len(args) != 3:
        raise InputError("Require InputGffFile, genomeSeqFastaFile and OutputFilename\n")
    if not os.path.exists(args[1]):
        raise InputError("Input file %s does not exist.\n" %args[1])
    
    #options.convertName = {}
    #if options.nameFile:
    #    f = open(options.nameFile, 'r')
    #    for line in f:
    #        items = line.strip().split()
    #        if len(items) < 2:
    #            continue
    #        options.convertName[items[0]] = items[1]
    #    f.close()

    genes = readGff( args[0] )
    header2seq = readFasta( args[1] )

    getGeneSeqs(genes, header2seq, args[2], options.qNameAttr)

if __name__ == '__main__':
    main()



















