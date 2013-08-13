#!/usr/bin/env python

'''
Input: Directory containing genes of each sample mapped to C.Ref. Format:
    indir/
        Sample1/
            file1.bed, ...
Output: Cluster all genes that were aligned together by Cactus. Out format:
<ID>\t<Sample1;Genes1>[\t<Samples2;Genes2>\t...]
'''

import os, sys, re, time, copy
from optparse import OptionParser

###################### Ojt classes #####################
class Bed():
    '''Bed record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 4: 
            raise BedFormatError("Bed format for this program requires a minimum of 4 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        self.chrom = items[0].split('.')[-1]
        self.chromStart = int(items[1]) #base 0
        self.chromEnd = int(items[2]) #exclusive
        self.name = items[3]
        self.bed12 = False

        if len(items) >= 12:
            self.bed12 = True
            self.score = int(items[4])
            self.strand = items[5]
            self.thickStart = int(items[6])
            self.thickEnd = int(items[7])
            self.itemRgb = items[8]
            self.blockCount = int(items[9])
            self.blockSizes = [ int(i) for i in items[10].split(',') ]
            self.blockStarts = [ int(i) for i in items[11].split(',') ]

    def __cmp__(self, other):
        if self.chrom != other.chrom:
            return cmp(self.chrom, other.chrom)
        elif self.strand != other.strand:
            return cmp(self.strand, other.strand)
        elif self.chromStart != other.chromStart:
            return cmp(self.chromStart, other.chromStart)
        else:
            return cmp(self.chromEnd, other.chromEnd)

class Region():
    '''A simple obj represents a region with chr, start, end
    '''
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.strand != other.strand:
            return cmp(self.strand, other.strand)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

class GeneFamily():
    '''Represent a gene family (containing a group of aligned genes)
    '''
    def __init__(self):
        self.beds = [] #list of beds. Each list represents one copy of the gene family mapped to CRef.
        self.sample2genes = {}


############### ERROR CLASSES #################
class BedFormatError(Exception):
    pass

############### FUNCTIONS #############
#============ misc. functions ================
def convertToAlnum(inputStr):
    items = re.split("[^a-zA-Z0-9]", inputStr)
    items = [item.capitalize() for item in items] 
    newStr = "".join(items)
    if not newStr.isalnum():
        print newStr
    return newStr

def getUnionList(list1, list2):
    unionset = set(list1) | set(list2)
    return list(unionset)

def getOtherList(mylist, mapToOther):
    #get the list of names in the other genome that are mapped to names in "mylist". 
    #Ignore names in mylist that have no mapping.
    otherlist = []
    for item in mylist:
        if item in mapToOther: #if gene of genome2 mapped to any gene of genome1
            otherlist.extend( mapToOther[item] )
    return otherlist

def getPc(count, total):
    if total == 0:
        return 0
    #return 100.0*count/total
    return 1.0*count/total

#============ Read input files ================
def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        items.append( line.strip() )
    f.close()
    #return list( set(items) )
    return items

def readBedFile(file):
    beds = {} #key = name (e.g geneName), val = list of beds for that gene
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        bed = Bed(line)
        if bed.name not in beds:
            beds[bed.name] = [bed]
        else:
            beds[bed.name].append(bed)
    f.close()
    return beds

def readParalogFile(file):
    fam2genes = {}
    gene2fam = {}
    #gene2paralogs = {} #key = gene, vals = list of paralogous genes
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        assert len(items) == 2
        genes = [ item.rstrip('|').split("|")[-1] for item in items[1].split(',') ]
        fam2genes[ items[0] ] = genes
        for g in genes:
            gene2fam[g] = items[0]
        #if len(genes) >= 2:
        #    for g in genes:
        #        gene2paralogs[ g ] = []
        #        for g2 in genes:
        #            if g != g2:
        #                gene2paralogs[g].append(g2)
    f.close()
    #return gene2paralogs
    return fam2genes, gene2fam

def readGeneList(file):
    genes = {} #key = gene, value = geneLenth
    f = open(file, 'r')
    for line in f:
        items = line.strip().split()
        if len(items) < 2:
            raise InputFormatError("genelist format has 2 fields: <geneNameOrID> <geneLength>. file %s, line %s only has %d fields\n" %(file, line, len(items)))
        gene = items[0]
        
        geneitems = gene.split('|')
        if len(geneitems) >=4:
            gene = geneitems[3]

        l = int(items[1])
        if gene not in genes:
            genes[gene] = l
        else:
            sys.stderr.write("Repetitive genes in genelist file %s, gene %s\n" %(file, gene))
    f.close()
    return genes

#============ comparing functions ================
def calcOverlap(bed1, bed2):
    if bed1.chrom != bed2.chrom or bed1.strand != bed2.strand or bed1.chromEnd <= bed2.chromStart or bed2.chromEnd <= bed1.chromStart:
        return 0.0, None
    intersectBed = copy.copy(bed2)
    blockSizes = []
    blockStarts = []

    for i1, start1 in enumerate( bed1.blockStarts ):
        start1 += bed1.chromStart
        end1 = start1 + bed1.blockSizes[i1]
        for i2, start2 in enumerate( bed2.blockStarts ):
            start2 += bed2.chromStart
            end2 = start2 + bed2.blockSizes[i2]
            if start1 < end2 and start2 < end1: #overlap
                start = max(start1, start2)
                end = min(end1, end2)
                blockStarts.append(start)
                blockSizes.append(end - start)
    
    overlap = sum(blockSizes)
    if overlap == 0:
        return 0.0, None

    intersectBed.chromStart = blockStarts[0]
    intersectBed.chromEnd = blockStarts[-1] + blockSizes[-1]
    intersectBed.thickStart = intersectBed.chromStart
    intersectBed.thickEnd = intersectBed.chromEnd
    intersectBed.blockCount = len( blockStarts )
    intersectBed.blockStarts = [s - intersectBed.chromStart for s in blockStarts]
    intersectBed.blockSizes = blockSizes
    return overlap, intersectBed

def checkHomologous( beds1, beds2, length1, coverage ):
    #calculate the number of shared bases between beds1 and beds2. If >= coverage, return the intersect beds
    overlap = 0.0
    beds = []
    for b1 in beds1:
        for b2 in beds2:
            o, b = calcOverlap(b1, b2)
            if o > 0:
                beds.append(b)
                overlap += o
    pc =getPc(overlap, length1)
    if pc >= coverage:#homologous, update beds2
        return beds
    else:
        return None

def addGeneToFams(fams, gene, sample, beds, genelen, coverage):
    for fam in fams:
        if famCheckHomologous(fam, beds, genelen, coverage): #found the orthologous family that gene belongs to
            fam.sample2genes[sample] = [gene]
            #HACK
            #if sample == "EscherichiaColi042Uid161985":
            #    print gene, beds
            #    print fam.sample2genes
            #    print fam.beds
            #END HACK
            return fam
    #Belong to a new gene family:
    fam = GeneFamily()
    fam.sample2genes[sample] = [gene]
    fam.beds = [beds]
    fams.append(fam)
    return fam

def famCheckHomologous(fam, beds, genelen, coverage):
    for i, currbeds in enumerate( fam.beds ):
        newbeds = checkHomologous(beds, currbeds, genelen, coverage)
        if newbeds:#found the aligned homolog
            fam.beds[i] = newbeds
            return True
    return False

def updateGeneFam(fam, gene, sample, beds, genelen, coverage):
    #add current gene to the gene family "fam" if >= coverage of its length is overlapped with fam beds
    fam.sample2genes[sample].append(gene)
    foundHomolog = famCheckHomologous(fam, beds, genelen, coverage)
    if not foundHomolog:
        fam.beds.append(beds)

def addSample(fams, gene2fam, sample, gene2beds, sampleGene2fam, sampleFam2genes, gene2len, coverage):
    #add a new sample and update the number of genes in the pan and core genomes
    addedSampleFams = 0
    for gene, beds in gene2beds.iteritems():
        assert gene not in gene2fam
        paralogs = []
        if gene in sampleGene2fam:
            paralogs = sampleFam2genes[ sampleGene2fam[gene] ]
        foundFam = False
        for p in paralogs:
            if p == gene:
                continue
            if p in gene2fam: #one of the current gene's paralogs has been visited - add this gene to that GeneFamily
                fam = gene2fam[p]
                updateGeneFam(fam, gene, sample, beds, gene2len[gene], coverage)
                gene2fam[gene] = fam
                foundFam = True
                break
        if not foundFam:
            fam = addGeneToFams(fams, gene, sample, beds, gene2len[gene], coverage)
            gene2fam[gene] = fam
            addedSampleFams += 1
    return addedSampleFams

def getOrthologGenes(indir, options):
    bedfiles = os.listdir(indir)
    if not options.samples:
        options.samples = bedfiles

    geneFams = []
    gene2fam = {}
    pans = []
    cores = []
    for i, sample in enumerate(options.samples):
        #get genelength:
        genelist = options.sample2genelist[sample]
        gene2len = readGeneList(genelist)

        #Read sample bed files
        sampledir = os.path.join(indir, sample)
        samplefiles = os.listdir( sampledir )
        sampleGene2beds = {}
        for samplefile in samplefiles:
            gene2beds = readBedFile(os.path.join(sampledir, samplefile))
            for g, beds in gene2beds.iteritems():
                sampleGene2beds[g] = beds
        #Paralog info:
        sampleFam2genes = {}
        sampleGene2fam = {}
        if sample in options.sample2gene2fam:
            sampleGene2fam = options.sample2gene2fam[sample]
            sampleFam2genes = options.sample2fam2genes[sample]

        addedFams = addSample(geneFams, gene2fam, sample, sampleGene2beds, sampleGene2fam, sampleFam2genes, gene2len, options.coverage)
        sampleSpecificFams = len(sampleFam2genes.keys()) - addedFams #gene family that is specific to the sample and therefore is missing from C.Ref
        
        #pan, core:
        pans.append( len(geneFams) + sampleSpecificFams )
        #print addedFams, sampleSpecificFams, pans[-1]
        if i == 0:
            cores.append( pans[0] )
        else:
            core = 0
            for fam in geneFams:
                if len(fam.sample2genes.keys()) == i + 1:
                    core += 1
            cores.append(core)
    return geneFams, pans, cores

def printOrthologGenes(outfile, geneFams):
    f = open(outfile, 'w')
    f.write("#ID\tSample1;genes1[\tSample2;genes2\t...]\n")
    for i, fam in enumerate(geneFams):
        f.write("%d" %i)
        for sample, genes in fam.sample2genes.iteritems():
            f.write("\t%s;%s" %(sample, ",".join(genes)))
        f.write("\n")
    f.close()

def printPanCore(outfile, samples, pans, cores):
    f = open(outfile, 'w')
    f.write("#AddedSample\tSampleCount\tPanGenome\tCoreGenome\tNewGeneFamilies\n")
    for i, sample in enumerate(samples):
        if i == 0:
            newfams = pans[0]
        else:
            newfams = pans[i] - pans[i-1]
        f.write("%s\t%d\t%d\t%d\t%d\n" %(sample, i + 1, pans[i], cores[i], newfams))
    f.close()

###################### MAIN ##################
def addOptions(parser):
    parser.add_option('-o', '--output', dest='outname', default='out', help='Output basename. Default=%default')
    parser.add_option('-c', '--coverage', dest='coverage', type='int', default=0.9, help='Minimum coverage for the "imperfect" aligned agreement category. Values: 0-1. Default=%default.')
    parser.add_option('--sampleOrder', dest='sampleOrder', help='The order in which the samples will be added in one at a time when calculating pan and core genome')
    parser.add_option('--paralog', dest='paralog', help='Directory containing sample paralog files, one file per sample. File format:<id>\t<comma,sep,list,of,paralogous,genes> ')
    parser.add_option('--genelist', dest='genelist', help='Directory containing <geneName>\t<geneLength> of each sample')

def checkOptions(parser, options, args):
    if len(args) < 1:
        parser.error("The input directory is required. None was given.")
    if not os.path.exists(args[0]) or not os.path.isdir(args[0]) :
        parser.error("Input directory %s does not exist or is not a directory.\n" %args[0])
    options.samples = None
    if options.sampleOrder:
        if not os.path.exists(options.sampleOrder):
            parser.error("sampleOrder file %s does not exist.\n" %(options.sampleOrder))
        else:
            options.samples = readList(options.sampleOrder)
    
    #options.sample2gene2paralogs = {}
    options.sample2fam2genes = {}
    options.sample2gene2fam = {}
    if options.paralog:
        files = os.listdir(options.paralog)
        for file in files:
            sample = convertToAlnum(file)
            #options.sample2gene2paralogs[sample] = readParalogFile( os.path.join(options.paralog, file) )
            options.sample2fam2genes[sample], options.sample2gene2fam[sample] = readParalogFile( os.path.join(options.paralog, file) )
            #print sample
            #print len(options.sample2fam2genes[sample])
    
    options.sample2genelist = {}
    if not options.genelist:
        parser.error("Genelist argument is required\n")
    else:
        files = os.listdir(options.genelist)
        for file in files:
            sample = convertToAlnum(file)
            options.sample2genelist[sample] = os.path.join(options.genelist, file)

def main():
    usage = "%prog <beddir>"
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, options, args)

    indir = args[0]
    geneFams, pans, cores = getOrthologGenes(indir, options)
    printOrthologGenes(options.outname, geneFams)
    printPanCore("%s-panNcore.txt" %options.outname, options.samples, pans, cores)

if __name__ == '__main__':
    main()

