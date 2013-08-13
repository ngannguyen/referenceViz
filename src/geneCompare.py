#!/usr/bin/env python

'''
Mon Jun 24 14:24:16 PDT 2013

Compare gene pairwise alignment with the "annotated" alignment.
Input: ( 0/ a query sample and a target sample)
         1/ List of annotated genes of the query sample. Format <genename>\t<genelength>
         2/ Directory contains annotated alignments of the query sample's genes with their homologs on the target sample.
            Format: indir1/
                        gene1.psl, gene2.psl, etc
            Where each psl entry represents alignment of the query gene to 1 homolog on the target sample. I.e if gene1
            has N homologs then file gene1.psl has N psl entries (lines)
         3/ Directory contains the gene alignments by an Aligner (e.g Cactus) to be compared with the annotated alignments
            Format: indir2/
                        gene1.psl, gene2.psl, etc
            Here, in case of rearrangement etc..., the alignment of 1 gene can span multiple psl entries. (I.e if gene1 
             maps to N places on the target, file gene1.psl can still have more than N psl entries in case of rearrangement/ bad alignments)
            Note: all the psl entries of each file are non-mergeable (see /hive/users/nknguyen/reconGit/mhc/bin/pslMerge.py)

Output: Define a query gene as mapping well to the target if the alignment mathces are >= 90% of the gene coverage and indels are of multiple of 3
        A/ Query sample's genes that map well to the target sample by the Aligner but do not have annotated homologs
            See why? Potential homologs? pseudogenes?
        
        B/ Query sample's genes with annotated homologs:
            1/ With 100% agreement between Aligner alignments & Annotated alignments
            2/ >= 90% agreement between Aligner alignments & Annotated alignments
            3/ Aligner maps the query genes well to the target sample, but to different places/locations compared with the Annotated alignments
               (unannotated multiple copies on the target??)
            4/ Aligner did not map the query genes well to the target (map partially or unmap).
'''

import os, sys, re, copy
from optparse import OptionParser

######################### Obj classes #####################
class Psl():
    '''Psl record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) != 21: 
            raise PslFormatError("Psl format requires 21 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        self.desc = line
        self.matches = int(items[0])
        self.misMatches = int(items[1])
        self.repMatches = int(items[2])
        self.nCount = int(items[3])
        self.qNumInsert = int(items[4])
        self.qBaseInsert = int(items[5]) #number of bases inserted in query
        self.tNumInsert = int(items[6]) #number of inserts in target
        self.tBaseInsert = int(items[7]) #number of bases inserted in target
        self.strand = items[8] #query strand
        self.qName = items[9] 
        self.qSize = int(items[10])
        self.qStart = int(items[11]) #base 0
        self.qEnd = int(items[12])
        self.tName = items[13]
        self.tSize = int(items[14])
        self.tStart = int(items[15])
        self.tEnd = int(items[16])
        self.blockCount = int(items[17])
        self.blockSizes = [int(s) for s in items[18].rstrip(',').split(',')]
        self.qStarts = [int(s) for s in items[19].rstrip(',').split(',')]
        self.tStarts = [int(s) for s in items[20].rstrip(',').split(',')]
        if len(self.blockSizes) != self.blockCount or len(self.qStarts) != self.blockCount or len(self.tStarts) != self.blockCount:
            raise PslFormatError("Psl format requires that the number of items in blockSizes, qStarts, tStarts is equal to blockCount. Line: %s\n" %line)

    def getDesc(self):
        blockSizes = ','.join([ str(s) for s in self.blockSizes])
        qStarts = ','.join([ str(s) for s in self.qStarts])
        tStarts = ','.join([ str(s) for s in self.tStarts])
        return "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s,\t%s,\t%s," %(self.matches, self.misMatches, self.repMatches, self.nCount, self.qNumInsert, self.qBaseInsert, self.tNumInsert, self.tBaseInsert, self.strand, self.qName, self.qSize, self.qStart, self.qEnd, self.tName, self.tSize, self.tStart, self.tEnd, self.blockCount, blockSizes, qStarts, tStarts)

    def reverse(self):
        #Reverse the strand of the query & target
        switchStrand = {'+': '-', '-': '+'}
        strand = ''
        for s in self.strand:
            if s in switchStrand:
                strand += switchStrand[s]
            else:
                strand += s
        self.strand = strand

        blockSizes = []
        qStarts = []
        tStarts = []
        for i in xrange(self.blockCount -1, -1, -1):
            blockSize = self.blockSizes[i]
            qEnd = self.qStarts[i] + blockSize
            rQstart = self.qSize - qEnd
            qStarts.append(rQstart)
            tEnd = self.tStarts[i] + blockSize
            rTstart = self.tSize - tEnd
            tStarts.append(rTstart)
            blockSizes.append(blockSize)
        self.blockSizes = blockSizes
        self.qStarts = qStarts
        self.tStarts = tStarts

######## ERROR CLASSES #######
class PslFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

######## COMPARE GENE ALIGNMENTS TO ANNOTATIONS #######
def getGoodAlignments(psls, length, minCoverage, indelMultipleOf3):
    #Return the psl that has total # of query bases = length & >= minCoverage matches & indel of multiple of 3
    assert length > 0
    goodPsls = []
    for psl in psls:
        #if indelMultipleOf3 and (psl.qBaseInsert % 3 != 0 or psl.tBaseInsert % 3 != 0): #indels not a multiple of 3
        #    continue
        matches = psl.matches + psl.repMatches
        aligned = matches + psl.misMatches
        #pcMatches = float(matches)/length 
        pcMatches = float(aligned)/length  #FLAG - may want to change this 
        if pcMatches < minCoverage: #not pass minimum alignment coverage
            continue
        queryBases = matches + psl.misMatches + psl.qBaseInsert
        if queryBases < length: #the psl entry does not represent the whole gene --> implies rearrangment since all the psls are un-mergeable
            continue
        if indelMultipleOf3 and aligned %3 != 0:#not a multiple of 3
            continue
        #Pass all the requirements:
        goodPsls.append(psl)
    return goodPsls

def unAnnotatedHomologs(annoGene2psls, gene2psls, gene2len, minCoverage):
    #Query sample's genes that map well to the target sample by the Aligner but do not have annotated homologs
    #See why? Potential homologs? pseudogenes?
    gene2psls_unAnnoGoodAlign = {}
    for gene, psls in gene2psls.iteritems():
        if gene not in annoGene2psls:
            assert gene in gene2len
            totalbases = gene2len[gene]
            #only care if the gene maps "well" to the target
            indelMul3 = True
            goodPsls = getGoodAlignments(psls, totalbases, minCoverage, indelMul3)
            if len(goodPsls) > 0:
                gene2psls_unAnnoGoodAlign[gene] = goodPsls
    return gene2psls_unAnnoGoodAlign

def printUnAnnoGoodAlign(gene2psls, outdir): 
    outfile = os.path.join(outdir, "unAnnotatedHomologs.txt")
    f = open(outfile, 'w')
    f.write("#Count: %d\n" %len(gene2psls))
    for gene, psls in gene2psls.iteritems():
        f.write("#Gene: %s\n" %gene)
        for psl in psls:
            f.write("%s\n" %psl.desc)
    f.close()

#========= comparing the homologs ========
def getAlignedPositions(psl):
    query2target = {} #key = position in query, value = position in target that aligns to key
    for i in xrange(psl.blockCount):
        qstart = psl.qStarts[i]
        tstart = psl.tStarts[i]
        size = psl.blockSizes[i]
        for j in xrange(size):
            qpos = qstart + j
            tpos = tstart + j
            query2target[qpos] = tpos
    return query2target

def getCommonAlignment(psl1, psl2):
    #return the number of aligned bases that are in both psl1 and psl2
    alignedBases = 0
    if psl1.qName != psl2.qName or psl1.tName != psl2.tName:
        return alignedBases
    aligned1 = getAlignedPositions(psl1)
    aligned2 = getAlignedPositions(psl2)
    for q1, t1 in aligned1.iteritems():
        if q1 in aligned2 and t1 == aligned2[q1]:
            alignedBases += 1
    return alignedBases

def compareAlignments(annoPsls, psls, minAgreement):
    psl2category = {} #categories: 'PA', 'IPA', 'D'
    for psl in psls: #for each good alignment
        alignedBases = psl.matches + psl.misMatches + psl.repMatches
        for annoPsl in annoPsls: #each annotated homolog
            annoAlignedBases = annoPsl.matches + annoPsl.misMatches + annoPsl.repMatches
            minAlignedBases = min([alignedBases, annoAlignedBases])
            if minAlignedBases == 0:
                continue

            numCommonAligned = getCommonAlignment(psl, annoPsl)
            pcAgreement = float(numCommonAligned)/minAlignedBases
            if pcAgreement >= minAgreement:
                if pcAgreement == 1.0:
                    psl2category[psl] = 'PA'
                else:
                    psl2category[psl] = 'IPA'
                break
        if psl not in psl2category: #did not agree with any of the annotated homologs
            psl2category[psl] = 'D'
    return psl2category

def compareHomologs(annoGene2psls, gene2psls, gene2len, minCoverage, minAgreement):
    #For query sample's genes with annotated homologs:
    #    1/ With 100% agreement between Aligner alignments & Annotated alignments
    #    2/ >= 90% agreement between Aligner alignments & Annotated alignments
    #    3/ Aligner maps the query genes well to the target sample, but to different places/locations compared with the Annotated alignments
    #       (unannotated multiple copies on the target??)
    #    4/ Aligner did not map the query genes well to the target (map partially or unmap).

    gene2psl2category = {} #categories: 'PA' (perfect agreement), 'IPA' (imperfect agreement), 'D' (disagreement), (if gene without psl --> 'U' (unmapped) category)
    for gene, psls in gene2psls.iteritems():
        if gene not in annoGene2psls:
            continue
        assert gene in gene2len
        genelen = gene2len[gene]
        indelMul3 = False
        goodPsls = getGoodAlignments(psls, genelen, minCoverage, indelMul3)
        annoPsls = annoGene2psls[gene]
        gene2psl2category[gene] = compareAlignments(annoPsls, goodPsls, minAgreement)
    return gene2psl2category

def getCategory2count(gene2psl2category):
    category2count = {'PA':0, 'IPA':0, 'D':0, 'U':0}
    for gene, psl2category in gene2psl2category.iteritems():
        if len(psl2category) == 0:
            category2count['U'] += 1
        else:
            cats = psl2category.values()
            if 'PA' in cats:
                category2count['PA'] += 1
            elif 'IPA' in cats:
                category2count['IPA'] += 1
            else:
                category2count['D'] += 1
    return category2count
            
def printCompareHomologs(gene2psl2category, outdir): 
    category2count = getCategory2count(gene2psl2category)
    numgenes = len(gene2psl2category)
    if numgenes == 0:
        return
    outfile = os.path.join(outdir, "compareHomologs.txt")
    f = open(outfile, 'w')
    f.write("#Genes with annotated homologs: %d\n" % numgenes)
    f.write("#Category\tCount\tPercentage\n")
    for c in ['PA', 'IPA', 'D', 'U']:
        count = category2count[c]
        pc = count*100.0/numgenes
        f.write("#%s\t%d\t%.3f\n" %(c, count, pc))

    f.write("#\n")
    for gene, psl2cat in gene2psl2category.iteritems():
        f.write("#Gene: %s\n" %gene)
        if len(psl2cat) == 0:
            f.write("#Yes_No\n")
        for psl, cat in psl2cat.iteritems():
            f.write("#%s\n" %cat)
            f.write("%s\n" %psl.desc)
    f.close()

def printSummaryStats(numTotal, numAnno, numAlign, numUnAnnoGood, gene2psl2category, outfile):
    category2count = getCategory2count(gene2psl2category)
    numUnAnno = numTotal - numAnno
    yes_yes = category2count['PA'] + category2count['IPA'] + category2count['D']
    numGoodAlign = numUnAnnoGood + yes_yes
    no_no = numUnAnno - numUnAnnoGood
    no_yes = numUnAnnoGood
    yes_no = category2count['U']

    stats = [numTotal, numAnno, numGoodAlign, no_no, no_yes, yes_no, yes_yes, category2count['PA'], category2count['IPA'], category2count['D']] 
    
    f = open(outfile, 'w')
    f.write("#Stats\tTotal\tAnnotated_Homologs\tCalled_Homologs\tNo_No\tNo_Yes\tYes_No\tYes_Yes\tPA\tIPA\tD\n")
    f.write("Count\t%s\n" % "\t".join([str(c) for c in stats]) )
    f.write("Count/Total\t%s\n" % "\t".join([ "%.3f" % (100.0*c/numTotal) for c in stats ]) )
    f.write("Count/Annotated_Homologs\tNA\t%s\n" % "\t".join(["%.3f" % (100.0*c/numAnno) for c in stats[1:]]) )
    f.write("Count/Called_Homologs\tNA\t%s\n" % "\t".join(["%.3f" % (100.0*c/numGoodAlign) for c in stats[1:]]) )
    f.write("Count/yes_yes\tNA\tNA\tNA\tNA\tNA\tNA\t%s\n" % "\t".join(["%.3f" % (100.0*c/yes_yes) for c in stats[6:]]) )
    f.close()

######## PROCESS INPUT FILES ########
def readPslFile(file):
    psls = [] 
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        psl = Psl(line)
        if len( psl.strand ) == 2 and psl.strand[1] == '-':
            psl.reverse()
        psls.append(psl)
    f.close()
    return psls

def readPslFiles(indir):
    gene2psls = {}
    for gene in os.listdir(indir):
        if gene.split('.')[-1] != 'psl':
            continue
        filename = os.path.join(indir, gene)
        psls = readPslFile(filename)
        gene = gene.rstrip('psl').rstrip('.')
        assert gene not in gene2psls
        gene2psls[gene] = psls
    return gene2psls

def readGeneList(file):
    gene2len = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        if len(items) < 2:
            raise InputError("Query sample gene list file has 2 fields with the format \
                              <gene>\\t<length>. File %s, line %s has only %d field(s).\n"\
                              %(file, line, len(items)))
        gene = items[0]
        length = int(items[1])
        if gene in gene2len:
            raise InputError("Redundant gene %s in file %s\n" %(gene, file))
        gene2len[gene] = length
    f.close()
    return gene2len

#######
def addOptions(parser):
    parser.add_option('-o', '--outdir', dest='outdir', help='Output directory')
    parser.add_option('--minCoverage', dest='minCoverage', type='float', default=0.9,\
                       help='Minimum portion of the query gene mapped to the target\
                       to be considered "well" mapped. Default = %default')
    parser.add_option('--minAgreement', dest='minAgreement', type='float', default=0.9,\
                       help='Minimum agreement between the Aligner and the annotated \
                       alignment. Default = %default')

def checkOptions(parser, args, options):
    if len(args) < 3:
        parser.error("Require at least 3 input arguments. Only got %d\n" %len(args))
    for i in xrange(len(args)):
        if not os.path.exists(args[i]):
            parser.error("%s does not exist\n" %args[i])
    if not options.outdir:
        options.outdir = os.getcwd()
    elif not os.path.exists(options.outdir):
        parser.error("Output directory %s does not exist\n" %options.outdir)

def main():
    usage = "%prog <queryGeneList> <alignments1> <alignments2> [options]"
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args() 
    checkOptions(parser, args, options)

    queryGeneFile = args[0]
    alignDir1 = args[1]
    alignDir2 = args[2]

    gene2len = readGeneList(queryGeneFile) #Query annotated genes
    gene2psls1 =  readPslFiles(alignDir1) #Annotated alignments
    gene2psls2 = readPslFiles(alignDir2) #Alignments to compare with

    gene2psls_unAnnoGoodAlign = unAnnotatedHomologs(gene2psls1, gene2psls2, gene2len, options.minCoverage)
    printUnAnnoGoodAlign(gene2psls_unAnnoGoodAlign, options.outdir)      

    gene2psl2category = compareHomologs(gene2psls1, gene2psls2, gene2len, options.minCoverage, options.minAgreement)
    printCompareHomologs(gene2psl2category, options.outdir)

    outfile = os.path.join(options.outdir, "summary.txt")
    printSummaryStats(len(gene2len), len(gene2psls1), len(gene2psls2), len(gene2psls_unAnnoGoodAlign), gene2psl2category, outfile)

if __name__ == '__main__':
    main()

