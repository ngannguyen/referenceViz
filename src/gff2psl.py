#!/usr/bin/env python

#Fri Feb  8 13:06:06 PST 2013
#nknguyen soe ucsc edu
#Convert a gff file to psl file
#Input: gff file
#Output: psl file

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

class Psl():
    '''Psl record
    '''
    def __init__(self):
        self.matches = 0 #number of bases that match that aren't repeats
        self.misMatches = 0 #number of bases that don't match
        self.repMatches = 0 #number of bases that match but are part of repeats
        self.nCount = 0 #number of N's bases
        self.qNumInsert = 0 #number of inserts in query
        self.qBaseInsert = 0 #number of bases inserted in query
        self.tNumInsert = 0 #number of inserts in target
        self.tBaseInsert = 0 #number of bases inserted in target
        self.strand = '.' #query strand
        self.qName = '' 
        self.qSize = 0 
        self.qStart = 0 #base 0
        self.qEnd = 0
        self.tName = ''
        self.tSize = 0
        self.tStart = 0
        self.tEnd = 0
        self.blockCount = 0
        self.blockSizes = []
        self.qStarts = []
        self.tStarts = []

    def getStr(self):
        blockSizes = ','.join([ str(s)for s in self.blockSizes])
        qStarts = ','.join([str(s) for s in self.qStarts])
        tStarts = ','.join([str(s) for s in self.tStarts])
        return "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s" \
               %(self.matches, self.misMatches, self.repMatches, self.nCount, 
                 self.qNumInsert, self.qBaseInsert, self.tNumInsert, self.tBaseInsert, self.strand, 
                 self.qName, self.qSize, self.qStart, self.qEnd, 
                 self.tName, self.tSize, self.tStart, self.tEnd, 
                 self.blockCount, blockSizes, qStarts, tStarts)


######## ERROR CLASSES #######
class GffFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

####### functions #########
def convertGff2psl(gffs, chr2size, qNameAttr, convertName):
    psl = Psl()
    
    if len(gffs) == 0:
        raise ValueError("convertGff2psl: no gff to convert\n")
        
    #Target information
    psl.tName = gffs[0].chr
    if psl.tName not in chr2size:
        raise InputFormatError
    psl.tSize = chr2size[psl.tName]
    if psl.tName in convertName:
        psl.tName = convertName[psl.tName]

    psl.tStart = gffs[0].start - 1 #gff is base 1 while psl is base 0
    psl.tEnd = gffs[-1].end #gff is inclusive while psl is exclusive
    tstrand = '+'

    psl.blockCount = len(gffs)
    
    #Query information
    qstrand = gffs[0].strand
    psl.strand = "%s%s" %(qstrand, tstrand)

    if qNameAttr not in gffs[0].attr:
        raise ValueError("qNameAttr %s is not in field 9 of gff %s\n" % (qNameAttr, gffs[0].desc))
    psl.qName = gffs[0].attr[qNameAttr]
    if psl.qName in convertName:
        psl.qName = convertName[psl.qName]
    
    for gff in gffs:#each exon
        if gff.strand != qstrand:
            raise ValueError("Strand of exons are inconsistent. GFF is probably unsorted. Sorted GFF is required.\n")
        blockSize = gff.end - gff.start + 1
        psl.blockSizes.append( blockSize )
        tStart = gff.start - 1 #convert from base 1 (gff) to base 0 (psl)
        psl.tStarts.append( tStart )
        if len(psl.tStarts) > 1:
            tPrevEnd = psl.tStarts[-2] + psl.blockSizes[-2]
            tBaseInsert = tStart - tPrevEnd
            if tBaseInsert > 0:
                psl.tNumInsert += 1
                psl.tBaseInsert += tBaseInsert

    psl.qStart = 0 
    psl.qSize = sum(psl.blockSizes)
    psl.qEnd = psl.qStart + psl.qSize
    psl.qStarts = [psl.qStart]
    for s in psl.blockSizes[:-1]:
        psl.qStarts.append(psl.qStarts[-1] + s)
    psl.matches = sum(psl.blockSizes)

    return psl

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

def convert2psl(genes, chr2size, outfile, qNameAttr, convertName):
    ofh = open(outfile, 'w')
    for gene in genes:
        gffs = []
        if len(gene.transcripts) == 0:
            psl = convertGff2psl([gene], chr2size, qNameAttr, convertName)
            ofh.write("%s\n" % psl.getStr())
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

                psl = convertGff2psl(gffs, chr2size, qNameAttr, convertName)
                ofh.write("%s\n" % psl.getStr())
    ofh.close()

def readChr2size(file):
    chr2size = {}
    f = open(file, 'r')
    for line in f:
        items = line.strip().split('\t')
        if len(items) < 2:
            raise InputFormatError("Wrong chrom2size format. Correct format is: <chr>\\t<size>\n")
        if items[0] in chr2size:
            raise InputFormatError("Repetitive chromosome in chr2size: %s\n" %items[0])
        chr2size[items[0]] = int(items[1])
    f.close()
    return chr2size

####### main ########
def main():
    usage = '%prog <chr2size> <inGffFile> <outPslFile>'
    parser = OptionParser(usage = usage)
    parser.add_option('--queryNameAttribute', dest='qNameAttr', default='ID', help='The attribute to extract (use as) query (gene) name from field 9 of the gff. Default = %default')
    parser.add_option('--convertSeqName', dest='nameFile', help='Optional. File that mapped sequence (chromosome) name in the gff to a different name in the output file. Format: <seqNameInGff> <newName>')

    options, args = parser.parse_args()
    if len(args) != 3:
        raise InputError("Require InputGffFile and OutputFilename\n")
    if not os.path.exists(args[1]):
        raise InputError("Input file %s does not exist.\n" %args[1])
    
    options.convertName = {}
    if options.nameFile:
        f = open(options.nameFile, 'r')
        for line in f:
            items = line.strip().split()
            if len(items) < 2:
                continue
            options.convertName[items[0]] = items[1]
        f.close()

    chr2size = readChr2size(args[0])
    genes = readGff( args[1] )
    convert2psl(genes, chr2size, args[2], options.qNameAttr, options.convertName)

if __name__ == '__main__':
    from referenceViz.src.gff2psl import *
    main()



















