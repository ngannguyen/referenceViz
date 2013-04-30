#!/usr/bin/env python2.6

#Mon Feb 18 15:19:49 PST 2013
#nknguyen soe ucsc edu
#Getting statistics of (common) genes mapped to a Reference sequence
#
#Input:
#   1/ Directory containing bed files, each file contains genes 
#      of each sample mapped to the Reference
#   2/ Directory containing genelist files, each file contains 
#      the annotated genelist of each sample
#   3/ File contains operon sets
#
#Output:
#   1/ For each sample:
#       how many annotated genes does the Ref contain
#       how many mapped perfectly to Ref (no indels)
#   2/ Draw: gene to number of samples for (annotated) and (mapped
#            perfectly to Ref)
#   3/ For the core genes: check to see if they align together
#   4/ Operons: 
#
#   5/ (Contiguity for genes): for each sample:
#      randomly pick any two genes, check to see if their order and
#      orientation is conserved on Ref.
#

import os, sys, re, time, random
from optparse import OptionParser
import libPlotting as libplot
import matplotlib.pyplot as pyplot

######################### Obj classes #####################
class Bed():
    '''Bed record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 4: 
            raise BedFormatError("Bed format for this program requires 4 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        #self.chr = items[0]
        self.chr = items[0].split('.')[1]
        self.start = int(items[1]) #base 0
        self.end = int(items[2]) #exclusive
        self.name = items[3]

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
class BedFormatError(Exception):
    pass

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

def getMapStatus(beds, totalbases):
    #totalbases = total number of bases of the annotated (sample) gene
    #Assuming the the beds are sorted in the order according to the gene coordinate
    #Check to see if the gene mapped perfectly to the reference (no gap, all blocks have the same order)
    #
    # [1s-->-1e][2s-->---2e][3s--->-3e]  or [3s--<-3e][2s--<---2e][1s---<-1e] : perfect
    # [1s-->-1e]---[2s-->---2e]---[3s--->-3e]  or [3s--<-3e]--[2s--<---2e]-------[1s---<-1e] : insertion
    # [1s-->-1e]-[2s-->2e]--[3s->-3e]  or [3s--<-3e][2s--<-2e][1s-<-1e] : deletion (original gene has bases not mapped to the ref
    # [2s-->-2e]---[1s-->---1e]---[3s--->-3e]  or [3s-->-3e][2s--<---2e][1s---<-1e] : non-linear rearrangement
    if len(beds) == 0:
        raise ValueError("trying to get status of an empty gene!\n")

    #Perfect when: del=False & rearrangement=False & insertion=False & folded=False
    status = {'deletion':0, 'rearrangement':False, 'insertion':0, 'folded':False}
    #Check the total number of bases mapped to Ref.
    mapped = sum( [bed.end - bed.start for bed in beds] )
    if totalbases < mapped:
        sys.stderr.write("Gene: %s, Total base: %d, mapped: %d\n" %(beds[0].name, totalbases, mapped))
        status['rearrangment'] = True
        status['Folded'] = True
    #assert totalbases >= mapped     #HACK! this should be un-commented
    
    if mapped < totalbases: 
        status['deletion'] = totalbases - mapped #number of bases that did not map to Ref.
    if len(beds) == 1:
        return status

    #Check if all the gene portions mapped to the same Ref. sequence (chromosome)
    for bed in beds:
        if bed.chr != beds[0].chr:
            status['rearrangement'] = True
            return status

    #all pieces mapped to same Ref. chromosome
    #Now check the ordering of the mapping:
    currBed = beds[0]
    currDirection = ''
    for i in xrange(1, len(beds)):
        bed = beds[i]
        if bed.start > currBed.start:
            direction = '+'
        else:
            direction = '-'
        if currDirection != '' and direction != currDirection: #change block ordering
            status['rearrangement'] = True
            return status
        currDirection = direction
        currBed = bed

    #Ordering is OK, check for overlapping
    forward = True
    if beds[0].start > beds[1].start:
        forward = False
    currBed = beds[0]
    totalgap = 0
    for i in xrange(1, len(beds)):
        bed = beds[i]
        if forward:
            gap = bed.start - currBed.end
        else:
            gap = currBed.start - bed.end
        if gap < 0: #overlap
            status['folded'] = True
            return status
        totalgap += gap
    status['insertion'] = totalgap

    return status

def writeCoverageHeader(f, fields):
    f.write("#Sample\tTotal")
    for field in fields:
        f.write("\t%s\t%%%s" %(field, field))
    f.write("\n")

def writeCoverageStats(f, fields, sample, stats, anno):
    total = stats['Total']
    f.write("%s\t%d" %(sample, total))
    for field in fields:
        if re.match("Anno", field):
            category = field.lstrip("Anno")
            f.write("\t%d\t%.2f" %( anno[category], getPc(anno[category], total) ))
        else:
            f.write("\t%d\t%.2f" %( stats[field], getPc(stats[field], total) ))
    f.write("\n")

def writeImperfectGeneBeds(f, gene, beds, category, sample):
    for bed in beds:
        f.write("%s\t%d\t%d\t%s\t%s\t%s\n" %(bed.chr, bed.start, bed.end, gene, category, sample))

def coverage(sample2fam2gene2beds, sample2fam2genes, gene2fam, gene2status, outbasename):
    '''
    For each sample:
        how many annotated genes does the Ref contain (partially and fully)? (mapped)
        how many annotated genes does the Ref contain (100% coverage) (fully mapped: deletion=0).
        how many mapped perfectly to Ref (perfect mapped, subset of fully mapped: deletion=0, insertion=0, rearrangement=False, folded=False)
        Print out the genes with deletion > 0
                            with insertion > 0
                            has rearrangement
                            is folded onto itself (duplications)
    '''
    sortedSample2genes = sorted( [(s, fam2genes) for s, fam2genes in sample2fam2genes.iteritems()], key=lambda item:len(item[1]) )
    sortedSamples = [ item[0] for item in sortedSample2genes ]

    famOutfile = "%s-coverage-fam.txt" %outbasename
    geneOutfile = "%s-coverage-gene.txt" %outbasename
    imperfectOutfile = "%s-imperfectGenes.txt" %outbasename #print genes that mapped imperfectly to C.Ref
    ffh = open(famOutfile, 'w')
    gfh = open(geneOutfile, 'w')
    ifh = open(imperfectOutfile, 'w')
    
    fields=['Mapped', 'FullyMapped', 'Deletion', 'Perfect', 'Insertion', 'AnnoInsertion', 'Rearrangement', 'AnnoRearrangement', 'Folded', 'AnnoFolded']
    writeCoverageHeader(ffh, fields)
    writeCoverageHeader(gfh, fields)
    ifh.write("#Chr\tStart\tEnd\tGene\tCategory\tSample\n")
    
    famAll = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0} #counts of geneFamilies
    famAllPc = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0} 
    geneAll = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0} #counts of genes
    geneAllPc = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
    
    famAllAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
    famAllPcAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
    geneAllAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
    geneAllPcAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}

    for sample in sortedSamples:
        fam2genes = sample2fam2genes[sample] #annotated (from the gff file)
        fam2gene2beds = sample2fam2gene2beds[sample] #genes mapped to C.Ref.
        famSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        geneSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        
        famAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        geneAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        
        famSample['Total'] = len(fam2genes)
        famSample['Mapped'] = len(fam2gene2beds)

        for fam, gene2beds in fam2gene2beds.iteritems():#each family with gene(s) mapped to C.Ref.
            assert fam in fam2genes
            total = len(fam2genes[fam])
            mapped = len(gene2beds.keys())
            currAnnoCount = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
            currCount = {'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
            
            for gene, beds in gene2beds.iteritems():
                #if gene not in gene2fam or gene2fam[gene] not in fam2genes or gene in gene2status:
                if gene not in gene2fam or gene2fam[gene] not in fam2genes:
                    sys.stderr.write("Sample %s, gene %s is not in gene2fam or the geneList.\n" %(sample, gene))
                    total -= 1
                    mapped -= 1
                    continue
                elif gene in gene2status:
                    annostatus = gene2status[gene]
                    geneAnno[ annostatus ] += 1
                    currAnnoCount[ annostatus ] += 1
                    continue
                     
                genelength = fam2genes[fam][gene]
                status = getMapStatus(beds, genelength)
                 
                if status['deletion'] == 0:
                    currCount['FullyMapped'] += 1
                    if not status['rearrangement'] and not status['folded'] and status['insertion'] == 0:
                        currCount['Perfect'] += 1
                else:
                    currCount['Deletion'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Deletion-%d' %status['deletion'], sample)

                if status['insertion'] > 0:
                    currCount['Insertion'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Insertion-%d' %status['insertion'], sample)
                if status['rearrangement']: 
                    currCount['Rearrangement'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Rearrangement', sample)
                if status['folded']: 
                    currCount['Folded'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Folded', sample)
            
            geneSample['Total'] += total
            geneSample['Mapped'] += mapped

            for category in currCount.keys():
                count = currCount[category]
                geneSample[category] += count
                if category in ['FullyMapped', 'Perfect']:
                    if count >= 1:
                        famSample[category] += 1
                elif count == mapped: #if all the mapped genes of the family have deletion/insertion/rearrangement/isFolded, then the family belongs to that category
                    famSample[category] += 1
            
            for k, v in currAnnoCount.iteritems():
                if v == len(gene2beds):
                    famAnno[k] += 1

        #Print coverage stats for the sample:
        writeCoverageStats(ffh, fields, sample, famSample, famAnno)
        writeCoverageStats(gfh, fields, sample, geneSample, geneAnno)
        
        #Cumulative stats
        for k, v in famSample.iteritems():
            famAll[k] += v
            geneAll[k] += geneSample[k]
            if k != 'Total':
                famAllPc[k] += getPc(v, famSample['Total'])
                geneAllPc[k] += getPc(geneSample[k], geneSample['Total'])

        for k, v in famAnno.iteritems():
            famAllAnno[k] += v
            geneAllAnno[k] += geneAnno[k]
            famAllPcAnno[k] += getPc(v, famSample['Total'])
            geneAllPcAnno[k] += getPc(geneAnno[k], geneSample['Total'])
    
    #Average
    n = len(sortedSamples)
    if n > 0:
        ffh.write("Average\t%.2f" %( famAll['Total']/n ))
        for field in fields:
            if re.match("Anno", field):
                field = field.lstrip("Anno")
                ffh.write( "\t%.2f\t%.2f" %(famAllAnno[field]/n, famAllPcAnno[field]/n) )
                gfh.write( "\t%.2f\t%.2f" %(geneAllAnno[field]/n, geneAllPcAnno[field]/n) )
            else:
                ffh.write( "\t%.2f\t%.2f" %(famAll[field]/n, famAllPc[field]/n) )
                gfh.write( "\t%.2f\t%.2f" %(geneAll[field]/n, geneAllPc[field]/n) )
        ffh.write("\n")

    ffh.close()
    gfh.close()
    ifh.close()

def getNumsam2numgene(gene2sams):
    numsam2numgene = {}
    for gene, sams in gene2sams.iteritems():
        numsam = len(sams)
        if numsam in numsam2numgene:
            numsam2numgene[numsam] += 1
        else:
            numsam2numgene[numsam] = 1
    return numsam2numgene

def checkAlign_psl(gene2perfectPsls):
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
        #else:
            #sys.stdout.write("\nPerfect but not allAligned: gene %s\n" %gene)
            #for p in psls:
            #    sys.stdout.write("%s\n" %p.desc)
    return numsam2numgene

def checkAlign(id2perfectBedLists):
    numsam2numgene = {}
    for fam, bedLists in id2perfectBedLists.iteritems():
        allAlign = True
        start = bedLists[0][0].start
        end = bedLists[0][-1].end
        size = sum([bed.end - bed.start  for bed in bedLists[0]])

        for beds in bedLists:
            currsize = sum([bed.end - bed.start for bed in beds])
            currstart = beds[0].start
            currend = beds[-1].end
            if (currsize == size and (currstart != start or currend != end)) or \
               (currsize < size and (currstart < start or currend > end)) or \
               (currsize > size and (currstart > start or currend < end)):
               allAlign = False
               break
        if allAlign:
            numsam = len(bedLists)
            if numsam in numsam2numgene:
                numsam2numgene[numsam] += 1
            else:
                numsam2numgene[numsam] = 1
        #else:
        #    sys.stdout.write("\nPerfect but not allAligned: gene family with id %s\n" %fam)
    return numsam2numgene

def drawSharedGenes(xdata, cat2ydata, outbasename, options):
    options.out = outbasename
    fig, pdf = libplot.initImage(11.2, 10.0, options)
    axes = fig.add_axes( [0.12, 0.18, 0.85, 0.75] )
    
    axes.set_title("Shared gene families")
    lines = []
    linenames = sorted( cat2ydata.keys() )
    colors = libplot.getColors6()
    markers = ['^', 'o', 'd', 'p', 'v', '*']
    for i, cat in enumerate(linenames):
        ydata = cat2ydata[cat]
        l = axes.plot(xdata, ydata, color=colors[i], marker=markers[i], markeredgecolor=colors[i], linestyle='-')
        lines.append(l)
    axes.set_xlabel( 'Number of samples' )
    axes.set_ylabel( 'Number of gene families' )

    axes.yaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    axes.xaxis.grid(b=True, color="#CCCCCC", linestyle='-', linewidth=0.005)
    
    legend = pyplot.legend( lines, linenames, numpoints=1, loc='best' )
    legend.__drawFrame = False

    libplot.writeImage( fig, pdf, options)
    

def sharedGenes(sample2fam2gene2beds, sample2fam2genes, gene2fam, outbasename, options):
    '''
    Count/draw: gene to number of samples for (annotated) and (mapped perfectly to Ref)
    '''
    fam2samAnno = {} 
    for sample, fams in sample2fam2genes.iteritems():
        for fam in fams:
            if fam not in sample2fam2gene2beds[sample]:
                sys.stderr.write("Sample %s: gene family with ID %s is in genelist but not in beds\n" %(sample, fam))
                continue

            if fam not in fam2samAnno:
                fam2samAnno[fam] = [sample]
            elif sample not in fam2samAnno[fam]:
                fam2samAnno[fam].append(sample)
    s2gAnno = getNumsam2numgene(fam2samAnno)

    famgene2samRef = {}
    famgene2samRefFull = {}
    famgene2samRefPerfect = {}
    famgene2perfectGenes = {}

    for sample, fam2gene2beds in sample2fam2gene2beds.iteritems():
        for fam, gene2beds in fam2gene2beds.iteritems():
            if fam not in sample2fam2genes[sample]:
                sys.stderr.write("Sample %s: gene family %s is in beds but not in genelist\n" %(sample, fam))
                continue
            
            if fam not in famgene2samRef:
                famgene2samRef[fam] = [sample]
            elif sample not in famgene2samRef[fam]:
                famgene2samRef[fam].append(sample)
            
            #check for full & perfect:
            isfull = False
            perfectGene = None
                
            for gene, beds in gene2beds.iteritems():
                genelength = sample2fam2genes[sample][fam][gene]
                status = getMapStatus(beds, genelength) ## {'deletion':int, 'rearrangement':False, 'insertion':int, 'folded':False}
                if status['deletion'] == 0:
                    isfull = True
                    if not status['rearrangement'] and not status['folded'] and status['insertion'] == 0:
                        perfectGene = beds
                        break
                if perfectGene:
                    break
            if isfull:
                if fam not in famgene2samRefFull:
                    famgene2samRefFull[fam] = [sample]
                elif sample not in famgene2samRefFull[fam]:
                    famgene2samRefFull[fam].append(sample)
            if perfectGene:
                if fam not in famgene2samRefPerfect:
                    famgene2samRefPerfect[fam] = [sample]
                    famgene2perfectGenes[fam] = [perfectGene]
                elif sample not in famgene2samRefPerfect[fam]:
                    famgene2samRefPerfect[fam].append(sample)
                    famgene2perfectGenes[fam].append(perfectGene)
    
    s2gRef = getNumsam2numgene(famgene2samRef)
    s2gRefFull = getNumsam2numgene(famgene2samRefFull)
    s2gRefPerfect = getNumsam2numgene(famgene2samRefPerfect)
    s2gRefPerfectAligned = checkAlign(famgene2perfectGenes)

    #Print to output file:
    outfile = "%s-gene2numsam" %outbasename
    f = open(outfile, 'w')
    f.write("#Number of gene families\tNumber of samples, annotated\tC.Ref\tC.Ref, fully mapped\tC.Ref, perfectly mapped\tC.Ref, perfectly mapped and all aligned\n")
    
    cumulativeAnno = 0
    cumulativeRef = 0
    cumulativeRefFull = 0
    cumulativeRefPerfect = 0
    cumulativeRefPerfectAligned = 0
    assert len(sample2fam2genes) == len(sample2fam2gene2beds)
    xdata = xrange(len(sample2fam2genes), 0, -1)
    cat2ydata = {'Annotated':[], 'CRef':[], 'CRefFull':[], 'CRefPerfect':[], 'CRefPerfectAligned':[]}

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
        cat2ydata['Annotated'].append(cumulativeAnno)
        cat2ydata['CRef'].append(cumulativeRef)
        cat2ydata['CRefFull'].append(cumulativeRefFull)
        cat2ydata['CRefPerfect'].append(cumulativeRefPerfect)
        cat2ydata['CRefPerfectAligned'].append(cumulativeRefPerfectAligned)

    f.close()

    #Draw scatter plot
    drawSharedGenes(xdata, cat2ydata, outbasename, options)

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

def getPosByCoverage(chr2pos2cov, cov):
    chr2pos = {}
    for chr, pos2cov in chr2pos2cov.iteritems():
        chr2pos[chr] = {}
        for pos, c in pos2cov.iteritems():
            if c == cov:
                chr2pos[chr][pos] = 1
    return chr2pos

def getNumsam2codingBases(sample2beds, wiggleFile):
    numsam2bases = {}
    chr2pos2cov = readWiggleFile(wiggleFile)

    visitedChr2pos = {}
    for sample, gene2beds in sample2beds.iteritems():
        for gene, beds in gene2beds.iteritems():
            for bed in beds:
                if bed.chr not in chr2pos2cov:
                    print "bedchr %s not in wiggle" % bed.chr
                    continue
                elif bed.chr not in visitedChr2pos:
                    visitedChr2pos[bed.chr] = {}

                for pos in xrange(bed.start, bed.end):
                    if pos in chr2pos2cov[bed.chr] and pos not in visitedChr2pos[bed.chr]:
                        visitedChr2pos[bed.chr][pos] = 1
                        numsam = chr2pos2cov[bed.chr][pos]
                        if numsam not in numsam2bases:
                            numsam2bases[numsam] = 1
                        else:
                            numsam2bases[numsam] += 1
    return numsam2bases

def getCoreBases2(sample2beds, numsam, wiggleFile):
    chr2pos2cov = readWiggleFile(wiggleFile)
    chr2pos = getPosByCoverage(chr2pos2cov, numsam)

    coreBases = 0
    visitedChr2pos = {}
    for sample, gene2beds in sample2beds.iteritems():
        for gene, beds in gene2beds.iteritems():
            for bed in beds:
                if bed.chr not in chr2pos:
                    print "bedchr %s not in chr2pos" % bed.chr
                    continue
                elif bed.chr not in visitedChr2pos:
                    visitedChr2pos[bed.chr] = {}

                for pos in xrange(bed.start, bed.end):
                    if pos in chr2pos[bed.chr] and pos not in visitedChr2pos[bed.chr]:
                        visitedChr2pos[bed.chr][pos] = 1
                        coreBases += 1
    return coreBases

def getCoreBases(sample2beds, numsam):
    samples = sample2beds.keys()
    sam2index = {}
    for i, sample in enumerate(samples):
        sam2index[sample] = i

    pos2samples = {} #key = position on C.Ref, val = list of sample indices
    for sample, gene2beds in sample2beds.iteritems():
        index = sam2index[sample]
        for gene, beds in gene2beds.iteritems():
            for bed in beds:
                for pos in xrange(bed.start, bed.end):
                    posStr = "%s.%d" % (bed.chr, pos)
                    if posStr not in pos2samples:
                        pos2samples[posStr] = [index]
                    elif index not in pos2samples[posStr]:
                        pos2samples[posStr].append(index)
    coreBases = 0
    for pos, indices in pos2samples.iteritems():
        if len(indices) == numsam:
            coreBases += 1

    #DEBUG
    sys.stdout.write( "Total positions: %d\n" % len(pos2samples.keys()) )
    numsam2numpos = {}
    for pos, indices in pos2samples.iteritems():
        ns = len(indices)
        if ns in numsam2numpos:
            numsam2numpos[ns] += 1
        else:
            numsam2numpos[ns] = 1
    sys.stdout.write("Numsam\tNumPos\n")
    for ns in sorted(numsam2numpos.keys()):
        sys.stdout.write( "%d\t%d\n" %(ns, numsam2numpos[ns]) )
    #END DEBUG
    
    return coreBases

######## PROCESS INPUT FILES ########
def convertToAlnum(inputStr):
    items = re.split("[^a-zA-Z0-9]", inputStr)
    items = [item.capitalize() for item in items] 
    newStr = "".join(items)
    assert newStr.isalnum()
    return newStr

def readBedFile(file):
    beds = {} #key = geneName (or geneFamilyId), val = list of beds for that gene
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

def readBedFiles(indir):
    sample2beds = {}
    for sample in os.listdir(indir):
        samplepath = os.path.join(indir, sample)
        if os.path.isdir(samplepath):
            files = [os.path.join(samplepath, file) for file in os.listdir(samplepath)]
        else:
            files = [samplepath]
        sample2beds[sample] = {}
        for file in files:
            beds = readBedFile(file)
            for name, bedlist in beds.iteritems():
                if name in sample2beds[sample]:
                    sample2beds[sample][name].extend(bedlist)
                else:
                    sample2beds[sample][name] = bedlist
    return sample2beds

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
        genenames = items[1].split(',')
        genes = []
        for g in genenames:
            #gi|386604534|ref|YP_006110834.1|
            geneitems = g.split('|')
            if len(geneitems) >= 4:
                g = geneitems[3]
            genes.append(g)
        for g in genes:
            gene2fam[g] = famid
    f.close()
    return gene2fam

def readWiggleFile(file):
    #fixedStep chrom=refChr1 start=1 step=1
    #57
    #60
    
    chr2pos2cov = {}
    f = open(file, 'r')
    chr = 'chr'
    start = 0 
    step = 0

    for line in f:
        line = line.strip()
        items = line.split()
        if len(items) == 4: #new chromosome
            chr = items[1].split("=")[1]
            start = int( items[2].split("=")[1] ) - 1 # -1 is to convert from base 1 to base 0
            step = int( items[3].split("=")[1] )
            pos = -1
            assert chr not in chr2pos2cov
            chr2pos2cov[chr] = {}
        else:
            if pos == -1:
                pos = start
            else:
                pos += step
            chr2pos2cov[chr][pos] = int(line)
            
    f.close()
    return chr2pos2cov

def readOperonFile(file):
    #First line contains specie name
    #
    operons = {} #key = operon name, val = list of genenames
    f = open(file, 'r')
    species = f.readline().strip()
    for line in f:
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        assert len(items) >= 3
        name = items[0]
        genes = items[1]
        strand = items[2]
        operons[name] = (genes, strand)

    f.close()
    return species, operons

def filterAbberantGenes(sample2genes, g2status):
    newsample2genes = {}
    for sample, gene2list in sample2genes.iteritems():
        genes = {}
        for gene, list in gene2list.iteritems():
            if gene not in g2status:
                genes[gene] = list
        newsample2genes[sample] = genes
    return newsample2genes

def readGeneStatus(file):
    gene2status = {}
    f = open(file, 'r')
    for line in f:
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.strip().split("\t")
        assert len(items) >= 2
        #gene2status[items[0]] = items[1]
        gene2status[items[0]] = items[1].split('-')[0]
    f.close()
    return gene2status

def getSample2fam2genes(sample2genes, gene2fam):
    sam2fam2genes = {}
    for sample, genes in sample2genes.iteritems():
        sam2fam2genes[sample] = {}
        for gene, length in genes.iteritems():
            if gene not in gene2fam:
                #raise ValueError("gene %s is not found in the geneClusters file.\n" %gene)
                sys.stderr.write("gene %s is not found in the geneClusters file. Sample %s\n" %(gene, sample))
                continue
            fam = gene2fam[gene]
            if fam not in sam2fam2genes[sample]:
                sam2fam2genes[sample][fam] = {gene:length}
            else:
                sam2fam2genes[sample][fam][gene] = length
    return sam2fam2genes

def getSample2fam2gene2beds(sample2beds, gene2fam):
    sam2fam2gene2beds = {}
    for sample, gene2beds in sample2beds.iteritems():
        fam2gene2beds = {}
        for gene, beds in gene2beds.iteritems():
            if gene not in gene2fam:
                #raise ValueError("gene %s is not found in the geneClusters file. Sample %s\n" %(gene, sample))
                sys.stderr.write("gene %s is not found in the geneClusters file. Sample %s\n" %(gene, sample))
                continue
            fam = gene2fam[gene]
            if fam not in fam2gene2beds:
                fam2gene2beds[fam] = {gene: beds}
            else:
                fam2gene2beds[fam][gene] = beds
        sam2fam2gene2beds[sample] = fam2gene2beds
    return sam2fam2gene2beds

def main():
    #bed directory: bed files, each containing genes of each sample mapped to C.Ref.
    #gene list directory: genelist files, each containing genes of each sample. Format: each line contains a gene/protein name
    #geneClusters: file clusters genes into family. Format: <GeneFamilyId>\t<Comma,sep,list,of,genes>
    
    usage = "%prog <bed directory> <gene list directory> <output basename> <geneClusters> [<gene2sample psl directory>]"
    #geneClusters has the format: <clusterID>\t<comma,separated,list,of,genes,within,this,cluster>
    parser = OptionParser(usage = usage)
    parser.add_option('-n', '--numsamplings', dest='numsamplings', type='int', default=100000, help='Number of samplings to performed for the contiguity stats. Default=%default')
    parser.add_option('-s', '--geneStatus', dest='gstatus', help='Specify the list of annotated genes that were folded/rearranged')
    parser.add_option('-w', '--wiggle', dest='wiggle', help='Wiggle file showing the coverage of each position')
    parser.add_option('-o', '--operons', dest='operonFile', help='File containing operon set. Default=%default' )

    libplot.initOptions(parser)
    options, args = parser.parse_args()
    libplot.checkOptions(options, parser)

    if len(args) < 4:
        parser.error("Required at least 4 inputs\n")

    g2status = {}
    if options.gstatus:
        g2status = readGeneStatus(options.gstatus)
    
    sample2beds = readBedFiles(args[0]) #genes mapped to CRef
    sample2genes = readGeneLists(args[1]) #annotated genes
    
    #number of coding bases share by all samples
    #numsam = len( sample2beds.keys() )
    #if options.wiggle:
    #    #coreBases = getCoreBases2(sample2beds, numsam, options.wiggle)
    #    numsam2bases = getNumsam2codingBases(sample2beds, options.wiggle)
    #    for numsam in sorted( numsam2bases.keys() ):
    #        print "%d\t%d" %(numsam, numsam2bases[numsam])
    #else:
    #    coreBases = getCoreBases(sample2beds, numsam)
    #    print "Coding core bases (shared by coding bases of %d samples):\n\t%d" %(numsam, coreBases)
    
    #if len(g2status) > 0:
    #    sample2beds = filterAbberantGenes(sample2beds, g2status)
    #    sample2genes = filterAbberantGenes(sample2genes, g2status)
    
    #cluster genes into gene families:
    gene2fam = readGeneClusters(args[3])
    sam2fam2gene2beds = getSample2fam2gene2beds(sample2beds, gene2fam) 
    sam2fam2genes = getSample2fam2genes(sample2genes,gene2fam)
    
    coverage(sam2fam2gene2beds, sam2fam2genes, gene2fam, g2status, args[2])
    #sharedGenes(sam2fam2gene2beds, sam2fam2genes, gene2fam, args[2], options) 
    
    if options.operonFile:
        sample, operons = readOperonFile(options.operonFile)
        checkOperons(operons, sam2fam2gene2beds[sample], sam2fam2genes[sample], gene2fam, g2status, args[2])

    #if len(args) == 4:
    #    sample2psls_sample = readPslFiles(args[3])
    #    contiguity(sample2psls_sample, sample2psls, args[2], options.numsamplings)

if __name__ == "__main__":
    from referenceViz.src.geneStats import *
    main()






