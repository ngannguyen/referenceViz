#!/usr/bin/env python

#Mon Feb 18 15:19:49 PST 2013
#nknguyen soe ucsc edu
#Getting statistics of (common) genes mapped to a Reference sequence
#
#Input:
#   1/ Directory containing psl files, each psl contains genes 
#      of each sample mapped to the Reference
#   2/ Directory containing genelist files, each file contains 
#      the annotated genelist of each sample
#
#Output:
#   1/ For each sample:
#       how many annotated genes does the Ref contain
#       how many mapped perfectly to Ref (no indels)
#   2/ Draw: gene to number of samples for (annotated) and (mapped
#            perfectly to Ref)
#   3/ For the core genes: check to see if they align together
#   4/ Operons
#
#   5/ (Contiguity for genes): for each sample:
#      randomly pick any two genes, check to see if their order and
#      orientation is conserved on Ref.
#

import os, sys, re, time, random
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

######## ERROR CLASSES #######
class PslFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

######## STATS ######
def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def coverage(sample2fam2gene2psls, sample2fam2genes, gene2fam, outbasename):
    '''
    For each sample:
        how many annotated genes does the Ref contain (partially and fully)? (mapped)
        how many annotated genes does the Ref contain (100% coverage) (fully mapped)
        how many mapped perfectly to Ref (no indels) (perfect mapped, subset of fully mapped)
    '''
    sortedSample2genes = sorted( [(s, fam2genes) for s, fam2genes in sample2fam2genes.iteritems()], key=lambda item:len(item[1]) )
    sortedSamples = [ item[0] for item in sortedSample2genes ]

    outfile = "%s-coverage.txt" %outbasename
    f = open(outfile, 'w')
    f.write("#Sample\tTotal\tMapped\t%Mapped\tFullyMapped\t%FullyMapped\tPerfectMapped\t%PerfectMapped\n")
    
    allTotal = 0.0
    allMapped = 0.0
    allPcMapped = 0.0
    allFullyMapped = 0.0
    allPcFullyMapped = 0.0
    allPerfect = 0.0
    allPcPerfect = 0.0

    for sample in sortedSamples:
        fam2genes = sample2fam2genes[sample]
        total = len(fam2genes)
        fam2gene2psls = sample2fam2gene2psls[sample]
        mapped = len(fam2gene2psls)

        fullyMapped = 0
        perfect = 0
        for fam, gene2psls in fam2gene2psls.iteritems():#each gene family, is fullymapped if >= 1 gene fullymapped, is perfect if >= 1 gene mapped perfectly
            isFullyMapped = False
            isPerfect = False
            for gene, psls in gene2psls.iteritems():
                if gene not in gene2fam or gene2fam[gene] not in fam2genes:
                    sys.stderr.write("Sample %s, gene %s is not in gene2fam or the geneList\n" %(sample, gene))
                    continue
                for psl in psls:
                    if psl.matches == psl.qSize:
                        isFullyMapped = True
                        if psl.qNumInsert == 0 and psl.tNumInsert == 0:
                            isPerfect = True
                            break
                if isPerfect:
                    break
            if isFullyMapped:
                fullyMapped += 1
            if isPerfect:
                perfect += 1
        f.write("%s\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n" %(sample, total, mapped, getPc(mapped, total), fullyMapped, getPc(fullyMapped, total), perfect, getPc(perfect, total) ))
            
        #Cumulative stats
        allTotal += total
        allMapped += mapped
        allPcMapped += getPc(mapped, total)
        allFullyMapped += fullyMapped
        allPcFullyMapped += getPc(fullyMapped, total)
        allPerfect += perfect
        allPcPerfect += getPc(perfect, total)

    #Average
    n = len(sortedSamples)
    if n > 0:
        f.write("Average\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" %(allTotal/n, allMapped/n, allPcMapped/n, allFullyMapped/n, allPcFullyMapped/n, allPerfect/n, allPcPerfect/n))

    f.close()

def getNumsam2numgene(gene2sams):
    numsam2numgene = {}
    for gene, sams in gene2sams.iteritems():
        numsam = len(sams)
        if numsam in numsam2numgene:
            numsam2numgene[numsam] += 1
        else:
            numsam2numgene[numsam] = 1
    return numsam2numgene

def checkAlign(gene2perfectPsls):
    numsam2numgene = {}
    for gene, psls in gene2perfectPsls.iteritems():
        allAlign = True
        start = psls[0].tStart
        end = psls[0].tEnd
        size = psls[0].qSize

        for psl in psls:
            if (psl.qSize == size and (psl.tStart != start or psl.tEnd != end)) or \
               (psl.qSize < size and (psl.tStart < start or psl.tEnd > end)) or\
               (psl.qSize > size and (psl.tStart > start or psl.tEnd < end)):
                allAlign = False
                break
                
                
        if allAlign:
            numsam = len(psls)
            if numsam in numsam2numgene:
                numsam2numgene[numsam] += 1
            else:
                numsam2numgene[numsam] = 1
        else:
            sys.stdout.write("\nPerfect but not allAligned: gene %s\n" %gene)
            for p in psls:
                sys.stdout.write("%s\n" %p.desc)
    return numsam2numgene

def sharedGenes(sample2fam2gene2psls, sample2fam2genes, gene2fam, outbasename):
    '''
    Count/draw: gene to number of samples for (annotated) and (mapped perfectly to Ref)
    '''
    fam2samAnno = {} 
    for sample, fams in sample2fam2genes.iteritems():
        for fam in fams:
            if fam not in sample2fam2gene2psls[sample]:
                sys.stderr.write("Sample %s: gene family with ID %s is in genelist but not in psls\n" %(sample, fam))
                continue

            if fam not in fam2samAnno:
                fam2samAnno[fam] = [sample]
            elif sample not in fam2samAnno[fam]:
                fam2samAnno[fam].append(sample)
    s2gAnno = getNumsam2numgene(fam2samAnno)

    famgene2samRef = {}
    famgene2samRefFull = {}
    famgene2samRefPerfect = {}
    famgene2perfectPsls = {}

    for sample, fam2gene2psls in sample2fam2gene2psls.iteritems():
        for fam, gene2psls in fam2gene2psls.iteritems():
            if fam not in sample2fam2genes[sample]:
                sys.stderr.write("Sample %s: gene family %s is in psls but not in genelist\n" %(sample, fam))
                continue
            
            if fam not in famgene2samRef:
                famgene2samRef[fam] = [sample]
            elif sample not in famgene2samRef[fam]:
                famgene2samRef[fam].append(sample)
            
            #check for full & perfect:
            isfull = False
            perfectPsl = None
            for gene, psls in gene2psls.iteritems():
                for psl in psls:
                    if psl.matches == psl.qSize:
                        isfull = True
                        if psl.qNumInsert == 0 and psl.tNumInsert == 0:
                            perfectPsl = psl
                            break
                if perfectPsl:
                    break
            if isfull:
                if fam not in famgene2samRefFull:
                    famgene2samRefFull[fam] = [sample]
                elif sample not in famgene2samRefFull[fam]:
                    famgene2samRefFull[fam].append(sample)
            if perfectPsl:
                if fam not in famgene2samRefPerfect:
                    famgene2samRefPerfect[fam] = [sample]
                    famgene2perfectPsls[fam] = [perfectPsl]
                elif sample not in gene2samRefPerfect[fam]:
                    famgene2samRefPerfect[fam].append(sample)
                    famgene2perfectPsls[fam].append(perfectPsl)
    
    s2gRef = getNumsam2numgene(famgene2samRef)
    s2gRefFull = getNumsam2numgene(famgene2samRefFull)
    s2gRefPerfect = getNumsam2numgene(famgene2samRefPerfect)
    s2gRefPerfectAligned = checkAlign(famgene2perfectPsls)

    #Print to output file:
    outfile = "%s-gene2numsam" %outbasename
    f = open(outfile, 'w')
    f.write("#Number of gene families\tNumber of samples, annotated\tC.Ref\tC.Ref, fully mapped\tC.Ref, perfectly mapped\tC.Ref, perfectly mapped and all aligned\n")
    
    cumulativeAnno = 0
    cumulativeRef = 0
    cumulativeRefFull = 0
    cumulativeRefPerfect = 0
    cumulativeRefPerfectAligned = 0
    assert len(sample2fam2genes) == len(sample2fam2gene2psls)
    for s in xrange(len(sample2fam2genes), 0, -1):
        if s in s2gAnno:
            cumulativeAnno += s2gAnno[s]
        if s in s2gRef:
            cumulativeRef += s2gRef[s]
        if s in s2gRefFull:
            cumulativeRefFull += s2gRefFull[s]
        if s in s2gRefPerfect:
            cumulativeRefPerfect += s2gRefPerfect[s]
        if s in s2gRefPerfectAligned:
            cumulativeRefPerfectAligned += s2gRefPerfectAligned[s]
        f.write( "%d\t%d\t%d\t%d\t%d\t%d\n" %(s, cumulativeAnno, cumulativeRef, cumulativeRefFull, cumulativeRefPerfect, cumulativeRefPerfectAligned) )

    f.close()

def sameChrom(psls1, psls2):
    for p1 in psls1:
        for p2 in psls2:
            if p1.tName == p2.tName:
                return True
    return False

def sampleContiguity(gene2psls1, gene2psls2, numsamplings):
    correct = 0
    for i in xrange(numsamplings):
        #randomly choose 2 genes:
        genes = random.sample(gene2psls1.keys() , 2)
        #keep sampling until find a gene pair that present in gene2psls2
        while genes[0] not in gene2psls2 or genes[1] not in gene2psls2 or not sameChrom(gene2psls1[genes[0]], gene2psls1[genes[1]]):
            genes = random.sample(gene2psls1.keys(), 2)
        
        psls11 = gene2psls1[genes[0]]
        psls12 = gene2psls1[genes[1]]

        psls21 = gene2psls2[genes[0]]
        psls22 = gene2psls2[genes[1]]
        
        reserve = False
        for p11 in psls11:
            for p12 in psls12:
                for p21 in psls21:
                    for p22 in psls22:
                        if p21.tName == p22.tName and p11.matches == p21.matches and p12.matches == p22.matches: #fully mapped to Cref
                            if ((p11.strand[0] == p21.strand[0] and p12.strand[0] == p22.strand[0]) and \
                                (p11.tStart - p12.tStart)*(p21.tStart - p22.tStart) > 0 ) or \
                               ((p11.strand[0] != p21.strand[0] and p12.strand[0] != p22.strand[0]) and \
                                (p11.tStart - p12.tStart)*(p21.tStart - p22.tStart) < 0 ):
                               reserve = True
                               break
                    if reserve:
                        break
                if reserve:
                    break
            if reserve:
                break
        if reserve:
            correct += 1
        else:
            sys.stdout.write("\nContiguity was broken for genes %s and %s:\n" %(genes[0], genes[1]))
            for psls in [psls11, psls12, psls21, psls22]:
                for p in psls:
                    sys.stdout.write("%s\n" %p.desc)

    return getPc(correct, numsamplings)

def contiguity(sample2psls1, sample2psls2, outbasename, numsamplings):
    '''
    Contiguity for genes: for each sample:
        randomly pick any two genes, check to see if their order and orientation on the sample (sample2psls1) is conserved on Ref (sample2psls2).
        repeat "numsamplings" times
    '''
    outfile = "%s-contiguity" %outbasename
    f = open(outfile, 'w')
    for sample, gene2psls1 in sample2psls1.iteritems():
        correct = sampleContiguity(gene2psls1, sample2psls2[sample], numsamplings) 
        f.write("%s\t%.2f\n" %(sample, correct))
    f.close()

######## PROCESS INPUT FILES ########
def convertToAlnum(inputStr):
    items = re.split("[^a-zA-Z0-9]", inputStr)
    items = [item.capitalize() for item in items] 
    newStr = "".join(items)
    assert newStr.isalnum()
    return newStr

def readPslFile(file):
    psls = {} #key = qName, val = Psl
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        psl = Psl(line)
        if psl.qName not in psls:
            psls[psl.qName] = [psl]
        else: #only take the longest transcripts
            currqSize = psls[psl.qName][-1].qSize
            if psl.qSize < currqSize:
                continue
            elif psl.qSize == currqSize:
                psls[psl.qName].append(psl)
            else:
                psls[psl.qName] = [psl]
    f.close()
    return psls

def readPslFiles(indir):
    sample2psls = {}
    for sample in os.listdir(indir):
        if sample.split('.')[-1] != 'psl':
            continue
        filename = os.path.join(indir, sample)
        psls = readPslFile(filename)
        sample = sample.rstrip('psl').rstrip('.')
        
        if re.search('_', sample): #HACK
            sample = convertToAlnum(sample) #HACK

        sample2psls[sample] = psls
    return sample2psls

def readGeneList(file):
    genes = []
    f = open(file, 'r')
    for line in f:
        gene = line.strip()
        if gene not in genes:
            genes.append(gene)
    f.close()
    return genes

def readGeneLists(indir):
    sample2genes = {}
    for sample in os.listdir(indir):
        filename = os.path.join(indir, sample)
        genes = readGeneList(filename)
        sample = convertToAlnum(sample) 
        sample2genes[sample] = genes
    return sample2genes

def readGeneClusters(file):
    gene2fam = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        if len(items) < 2:
            continue
        famid = items[0]
        genes = items[1].split(',')
        #gi|386604534|ref|YP_006110834.1|
        genes = [g.split('|')[4] for g in genes]
        for g in genes:
            gene2fam[g] = famid
    f.close()
    return gene2fam

def getSample2Fam2Genes(sample2genes,gene2fam):
    sam2fam2genes = {}
    for sample, genes in sample2genes.iteritems():
        sam2fam2genes[sample] = {}
        for gene in genes:
            if gene not in gene2fam:
                raise ValueError("gene %s is not found in the geneClusters file.\n" %gene)
            fam = gene2fam[gene]
            if fam not in sam2fam2genes[sample]:
                sam2fam2genes[sample][fam] = [gene]
            else:
                sam2fam2genes[sample][fam].append(gene)
    return sam2fam2genes

def getSample2Fam2Gene2Psls(sample2psls, gene2fam):
    sam2fam2gene2psls = {}
    for sample, gene2psls in sample2psls.iteritems():
        sam2fam2gene2psls[sample] = {}
        fam2gene2psls = {}
        for gene, psls in gene2psls.iteritems():
            if gene not in gene2fam:
                raise ValueError("gene %s is not found in the geneClusters file.\n" %gene)
            fam = gene2fam[gene]
            if fam not in fam2gene2psls:
                fam2gene2psls[fam] = {gene: psls}
            else:
                fam2gene2psls[fam][gene] = psls
    return sam2fam2gene2psls

def main():
    usage = "%prog <psl directory> <gene list directory> <output basename> <geneClusters> [<gene2sample psl directory>]"
    #geneClusters has the format: <clusterID>\t<comma,separated,list,of,genes,within,this,cluster>
    parser = OptionParser(usage = usage)
    parser.add_option('-n', '--numsamplings', dest='numsamplings', type='int', default=100000, help='Number of samplings to performed for the contiguity stats. Default=%default')
    options, args = parser.parse_args()

    if len(args) < 4:
        parser.error("Required at least 4 inputs\n")

    sample2psls = readPslFiles(args[0])
    sample2genes = readGeneLists(args[1])
    
    gene2fam = readGeneClusters(args[3])
    sam2fam2gene2psls = getSample2Fam2Gene2Psls(sample2psls, gene2fam) 
    sam2fam2genes = getSample2Fam2Genes(sample2genes,gene2fam)

    coverage(sam2fam2gene2psls, sam2fam2genes, gene2fam, args[2])
    sharedGenes(sam2fam2gene2psls, sam2fam2genes, gene2fam, args[2]) 
    
    #if len(args) == 4:
    #    sample2psls_sample = readPslFiles(args[3])
    #    contiguity(sample2psls_sample, sample2psls, args[2], options.numsamplings)

if __name__ == "__main__":
    from referenceViz.src.pslGeneStats import *
    main()






