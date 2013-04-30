#!/usr/bin/env python

#Sun Mar 17 22:16:56 PDT 2013
#nknguyen soe ucsc edu
#Input: gff
#Output: bacterial genes (should be one exon) with abnormal status (folded / rearrangment/ insertion)
#(insertion: when gene has >=2 exons, status = "insertion-totalBasesBtwExons")
#<genename>\t<status>
#

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

    #Return only genes with CDS(s):
    genesWithCDSs = []
    for g in genes:
        if len(g.transcripts) > 0:
            genesWithCDSs.append(g)

    f.close()
    return genesWithCDSs

def getGeneStatus(gffs):
    #Check for same chrom:
    for gff in gffs:
        if gff.chr != gffs[0].chr:
            return 'Rearrangement'
    #Check the ordering and overlapping: make sure exon0.start <= exon0.end < exon1.start <= exon1.end < exon2.start ...
    currGff = gffs[0]
    assert currGff.start <= currGff.end
    introns = 0
    folded = 0
    for i in xrange(1, len(gffs)):
        gff = gffs[i]
        if currGff.end >= gff.start:
            folded += currGff.end - gff.start +1
            #return 'Folded'
        else:
            assert gff.start <= gff.end
            introns += gff.start - currGff.end - 1
        currGff = gff
    
    if introns == 0 and folded == 0:
        return "OK"    
    elif folded > 0:
        return 'Folded-%d' %folded
    else:
        return "Insertion-%d" %introns

def getStatus(genes, outfile, nameAttr, sampleName):
    ofh = open(outfile, 'w')
    if sampleName:
        ofh.write("#%s\n" %sampleName)
    else:
        ofh.write("#%s\tStatus\n" %nameAttr)
    pos2count = {}

    for gene in genes:
        gffs = []
        if len(gene.transcripts) > 0:
            for tx in gene.transcripts:
                #status = "OK"
                if len(tx.exons) > 1:
                    status = getGeneStatus( tx.exons )
                    if nameAttr not in tx.exons[0].attr:
                        if nameAttr not in gene.attr:
                            raise ValueError("nameAttr %s is not in field 9 of gff %s\n" %(nameAttr, gene.desc))
                        tx.exons[0].attr[ nameAttr ] = gene.attr[ nameAttr ]

                    if status != "OK":
                        ofh.write("%s\t%s\n" % (tx.exons[0].attr[nameAttr], status) )
                
                #if the transcript status is ok, record the position of the gene
                #if status == "OK":
                for exon in tx.exons:
                    for p in xrange(exon.start, exon.end + 1):
                        if p not in pos2count:
                            pos2count[p] = 1
                        else:
                            sys.stdout.write("gene %s, exon %s, start: %d, end: %d\n" %(gene.attr[nameAttr], exon.attr[nameAttr], exon.start, exon.end))
                            pos2count[p] += 1
    sys.stdout.write( "%d\n" % len(pos2count.keys()) )
    sys.stdout.write( "Redundant %d\n" % sum(pos2count.values()) )
    ofh.close()

####### main ########
def main():
    usage = '%prog <gff> <outputFile>'
    parser = OptionParser(usage = usage)
    parser.add_option('-a', '--queryNameAttribute', dest='qNameAttr', default='ID', help='The attribute to extract (use as) query (gene) name from field 9 of the gff. Default = %default')
    parser.add_option('-s', '--sampleName', dest='sampleName', help='Sample name (optional). Default = %default')

    options, args = parser.parse_args()
    if len(args) != 2:
        raise InputError("Require InputGffFile, and OutputFilename\n")
    if not os.path.exists(args[0]):
        raise InputError("Input file %s does not exist.\n" %args[0])
    
    genes = readGff( args[0] )
    getStatus(genes, args[1], options.qNameAttr, options.sampleName)

if __name__ == '__main__':
    main()



















