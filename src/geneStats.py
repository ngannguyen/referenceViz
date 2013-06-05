#!/usr/bin/env python2.6

#EDIT:
#Tue May 14 13:42:02 PDT 2013
#Mon May  6 12:44:59 PDT 2013
#Mon Feb 18 15:19:49 PST 2013
#nknguyen soe ucsc edu
#Getting statistics of (common) genes mapped to a Reference sequence
#
#Input:
#   1/ Directory containing bed files, each file contains genes 
#      of each sample mapped (lift-overed) to the Reference
#   2/ Directory containing genelist files, each file contains 
#      the annotated genelist of each sample
#   3/ geneClusters: file showing groups of homologous genes. Format: <groupID>\t<gene1,gene2,...>
#   
#   Optional:
#   1/ File contains operon sets
#   2/ refGeneBeds: bed file of the annotated genes of the Reference 
#
#Output:
#   Analysis 1: how genes of each sample mapped to the Reference
#       For each sample:
#           how many annotated genes does the Reference contain
#           how many mapped perfectly to Ref (no indels)
#
#   Analysis 2: how are the alignments used to lift over genes of each sample to the Reference
#               compared with the annotated gene sets (gene set of each sample and gene set of the Reference)
#
#   (Others:
#   3/ Draw: gene to number of samples for (annotated) and (mapped
#            perfectly to Ref)
#   4/ For the core genes: check to see if they align together
#   5/ Operons: 
#
#   6/ (Contiguity for genes): for each sample:
#      randomly pick any two genes, check to see if their order and
#      orientation is conserved on Ref.
#   )

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
        self.chr = items[0].split('.')[-1]
        self.start = int(items[1]) #base 0
        self.end = int(items[2]) #exclusive
        self.name = items[3]

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

######## ERROR CLASSES #######
class BedFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

class Region():
    '''A simple obj represents a region with chr, start, end
    '''
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

######## ERROR CLASSES #######
class BedFormatError(Exception):
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
    #status = {'deletion':0, 'rearrangement':False, 'insertion':0, 'folded':False}
    #Perfect when: del=0 & rearrangement=False & insertion=0 & folded=False
    #
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

def getNonOverlapRegions(beds):
    chr2regs = {}
    beds = sorted(beds)
    #sys.stderr.write("getNonOverlapRegions... Done sorting input beds\n")

    for bed in beds:
        reg = Region(bed.chr, bed.start, bed.end)
        if reg.chr not in chr2regs: #new chromosome
            chr2regs[reg.chr] = [reg]
        else:
            hasoverlap = False
            for r in chr2regs[reg.chr]:
                if reg.start < r.end and r.start < reg.end: #overlap
                    r.start = min( r.start, reg.start )
                    r.end = max( r.end, reg.end )
                    hasoverlap = True
                    break
            if not hasoverlap: #non-overlap
                chr2regs[reg.chr].append(reg)
    return chr2regs

def calcOverlap(regs1, regs2):
    assert len(regs1) > 0 and len(regs2) > 0
    i1 = 0
    i2 = 0
    overlap = 0.0
    while i1 < len(regs1) and i2 < len(regs2):
        curr1 = regs1[i1]
        curr2 = regs2[i2]
        if curr1.start < curr2.end and curr2.start < curr1.end: #overlap
            overlap += min(curr1.end, curr2.end) - max(curr1.start, curr2.start)
            if curr1.end < curr2.end:
                i1 += 1
            else: #curr2.end < curr1.end:
                i2 += 1
        else: #non overlap
            if curr1.end <= curr2.start:
                i1 += 1
            else:
                i2 += 1
    return overlap

def compareAlignments( beds1, beds2, coverage ):
    #compare two alignments, beds have only chr, start & end info
    #sys.stderr.write("comparing alignments. beds1 %s, %d. beds2 %s, %d\n" %(beds1[0].name, len(beds1), beds2[0].name, len(beds2)))
    chr2regs1 = getNonOverlapRegions(beds1) 
    chr2regs2 = getNonOverlapRegions(beds2) 
    #DEBUG
    #l1 = 0
    #for c, rs in chr2regs1.iteritems():
    #    l1 += len(rs)
    #l2 = 0
    #for c, rs in chr2regs2.iteritems():
    #    l2 += len(rs)
    #sys.stderr.write("nonOverlapRegs1: %d; regs2: %d\n" %(l1, l2))
    #END DEBUG

    total1 = 0.0
    total2 = 0.0
    overlap = 0.0
    for chr, regs1 in chr2regs1.iteritems():
        #total1 += sum([reg.end - reg.start for reg in regs1])
        for r1 in regs1:
            total1 += r1.end - r1.start

        if chr not in chr2regs2:
            continue
        
        regs2 = chr2regs2[chr]
        #total2 += sum([reg.end - reg.start for reg in regs2])
        for r2 in regs2:
            total2 += r2.end - r2.start
        
        overlap += calcOverlap(regs1, regs2) 
   
    pc1 = getPc(overlap, total1)
    pc2 = getPc(overlap, total2)
    if pc1 == 100.0 and pc2 == 100.0: 
        return 1 #perfect aligned
    elif pc1 >= coverage and pc2 >= coverage:
        return 2 #imperfect aligned
    else:
        return 3 #other

def checkAligned( qbeds, tgene2beds, tgenes ):
    #qbeds:
    minCoverage = 90.0
    #categories = {1: 'totalAligned', 2: 'partiallyAligned', 3:'other'}
    categories = {1: 'A', 2: 'IA', 3:'O'}

    cats = []
    for tgene in tgenes:
        if tgene not in tgene2beds:
            sys.stderr.write("tgene %s is not in tgene2beds\n" %tgene)
            cats.append(3)
            continue
        tbeds = tgene2beds[tgene]
        cats.append( compareAlignments(qbeds, tbeds, minCoverage) )

    return categories[ min(cats) ]

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

def writeAlignerComparisonStats(afh, fields, sample, alignerAgreement, total):
    #all_alignerAgreement = {'YesP':0.0, 'YesIns':0.0, 'YesDel':0.0, 'YesFolded': 0.0, 'YesRearrange':0.0, \
    #                    'NoP':0.0, 'NoIns':0.0, 'NoDel':0.0, 'NoFolded':0.0, 'NoRearrange':0.0, 'YesNo':0.0 } #ref_fam2genes
    #aFields = sorted( all_alignerAgreement.keys() )
    #afh.write( "Sample\t%s\tYesYes\tNoYes\n" %( "\t".join(aFields) ) )
    yy = 0.0
    ny = 0.0
    afh.write("%s\t%d" % (sample, total))
    for k in fields:
        v = alignerAgreement[k]
        afh.write("\t%d" %v)
        if re.match("Yes", k) and k != "YesNo" and k != "YesPAligned":
            yy += v
        if re.match("No", k):
            ny += v
    afh.write("\t%d\t%d\n" %(yy, ny))

def writeImperfectGeneBeds(f, gene, beds, category, sample):
    for bed in beds:
        f.write("%s\t%d\t%d\t%s\t%s\t%s\n" %(bed.chr, bed.start, bed.end, gene, category, sample))

def coverage2(refname, ref_gene2beds, sample2fam2gene2beds, sample2fam2genes, gene2fam, gene2status, outbasename):
    '''
    For each sample:
        how many annotated genes does the Ref contain (partially and fully)? (mapped)
        how many annotated genes does the Ref contain (100% coverage) (fully mapped: deletion=0).
        how many mapped perfectly to Ref (perfect mapped, subset of fully mapped: deletion=0, insertion=0, rearrangement=False, folded=False)
        Print out the genes with deletion > 0
                            with insertion > 0
                            that have rearrangement
                            that are folded onto themselves (duplications)
    '''
    #sort samples by the number of gene families
    sortedSample2genes = sorted( [(s, fam2genes) for s, fam2genes in sample2fam2genes.iteritems()], key=lambda item:len(item[1]) )
    sortedSamples = [ item[0] for item in sortedSample2genes ]
    for sample in sortedSamples:
        if sample not in sample2fam2gene2beds or sample == refname:
            sortedSamples.remove(sample)

    #Annotated genes of the reference genome:
    ref_fam2genes = sample2fam2genes[refname] 

    famOutfile = "%s-coverage-fam.txt" %outbasename
    geneOutfile = "%s-coverage-gene.txt" %outbasename
    alignerOutfile = "%s-alignerComparisons.txt" %outbasename
    imperfectOutfile = "%s-imperfectGenes.txt" %outbasename #print genes that mapped imperfectly to C.Ref
    
    ffh = open(famOutfile, 'w')
    gfh = open(geneOutfile, 'w')
    afh = open(alignerOutfile, 'w')
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

    #all_alignerAgreement = {'YesP':0.0, 'YesPAligned':0.0, 'YesIns':0.0, 'YesDel':0.0, 'YesFolded': 0.0, 'YesRearrange':0.0, \
    #                    'NoP':0.0, 'NoIns':0.0, 'NoDel':0.0, 'NoFolded':0.0, 'NoRearrange':0.0, 'YesNo':0.0 } #ref_fam2genes
    all_alignerAgreement = {'YesPA':0.0, 'YesPIA':0.0, 'YesPO':0.0, \
                            'YesInsIA':0.0, 'YesInsO':0.0, 'YesDelIA':0.0, 'YesDelO':0.0, \
                            'YesFoldedIA':0.0, 'YesFoldedO':0.0, 'YesRearrangeIA':0.0, 'YesRearrangeO':0.0, \
                            'NoP':0.0, 'NoIns':0.0, 'NoDel':0.0, 'NoFolded':0.0, 'NoRearrange':0.0, 'YesNo':0.0 }#ref_fam2genes, A: 100% aligned, IA: >= XX% aligned, O: others
    aFields = sorted( all_alignerAgreement.keys() )
    afh.write( "Sample\tTotal\t%s\tYesYes\tNoYes\n" %( "\t".join(aFields) ) )

    for sample in sortedSamples:
        fam2genes = sample2fam2genes[sample] #annotated (from the gff file)
        fam2gene2beds = sample2fam2gene2beds[sample] #genes mapped to Ref.
        famSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        geneSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        
        famAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        geneAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        
        famSample['Total'] = len(fam2genes)
        famSample['Mapped'] = len(fam2gene2beds)

        #alignerAgreement = {'YesP':0.0, 'YesPAligned':0.0, 'YesIns':0.0, 'YesDel':0.0, 'YesFolded': 0.0, 'YesRearrange':0.0, \
        #                    'NoP':0.0, 'NoIns':0.0, 'NoDel':0.0, 'NoFolded':0.0, 'NoRearrange':0.0, 'YesNo':0.0 } #ref_fam2genes
        #ref_fam2genes, A: 100% aligned, IA: >= XX% aligned, O: others
        alignerAgreement = {'YesPA':0.0, 'YesPIA':0.0, 'YesPO':0.0, \
                            'YesInsIA':0.0, 'YesInsO':0.0, 'YesDelIA':0.0, 'YesDelO':0.0, \
                            'YesFoldedIA':0.0, 'YesFoldedO':0.0, 'YesRearrangeIA':0.0, 'YesRearrangeO':0.0, \
                            'NoP':0.0, 'NoIns':0.0, 'NoDel':0.0, 'NoFolded':0.0, 'NoRearrange':0.0, 'YesNo':0.0 }
        #number of genes of "sample" aligned to at least one gene in Ref by BLAT but NOT by cactus:
        blatonly = 0.0
        geneSampleTotal = 0.0
        for fam, genes in fam2genes.iteritems():
            #if len( fam2genes[fam] ) > 1:#SKIP DUPLICATED GENES
            #    continue
            if fam not in fam2gene2beds and fam in ref_fam2genes:
                blatonly += len(genes)
            geneSampleTotal += len(genes)
        alignerAgreement[ 'YesNo' ] = blatonly

        for fam, gene2beds in fam2gene2beds.iteritems():#each family with gene(s) mapped to Ref.
            assert fam in fam2genes
            #if len( fam2genes[fam] ) > 1:#SKIP DUPLICATED GENES
            #    continue

            #total = len(fam2genes[fam])
            mapped = len(gene2beds.keys())
            currAnnoCount = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
            currCount = {'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
            
            blat = "No"
            if fam in ref_fam2genes: #there are genes in Ref aligned with genes of current genome
                blat = "Yes"

            for gene, beds in gene2beds.iteritems():#each gene in the family
                cactus = 'NP'

                #if gene not in gene2fam or gene2fam[gene] not in fam2genes or gene in gene2status:
                if gene not in gene2fam or gene2fam[gene] not in fam2genes:
                    sys.stderr.write("Sample %s, gene %s is not in gene2fam or the geneList.\n" %(sample, gene))
                    #total -= 1
                    geneSampleTotal -= 1
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
                        cactus = 'P'
                else:
                    currCount['Deletion'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Deletion-%d' %status['deletion'], sample)
                    cactus = 'Del'

                if status['insertion'] > 0:
                    currCount['Insertion'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Insertion-%d' %status['insertion'], sample)
                    cactus = 'Ins'
                if status['rearrangement']: 
                    currCount['Rearrangement'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Rearrangement', sample)
                    cactus = 'Rearrange'
                if status['folded']: 
                    currCount['Folded'] += 1
                    writeImperfectGeneBeds(ifh, gene, beds, 'Folded', sample)
                    cactus = 'Folded'

                blatcactus = "%s%s" % (blat, cactus)
                #alignerAgreement[ blatcactus ] += 1.0
                #if blatcactus == "YesP" :
                if blat == "Yes": #(cactus != "No")
                    alignedStatus = checkAligned(beds, ref_gene2beds, ref_fam2genes[fam]) #A: 100% aligned, IA (imperfect aligned): >= XX% aligned, O: others
                    if alignedStatus == 'A' and cactus != 'P':#HACK
                        blatcactus += 'IA'
                    else:
                        blatcactus += alignedStatus
                alignerAgreement[ blatcactus ] += 1.0
            
            #geneSample['Total'] += total
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
        geneSample['Total'] = geneSampleTotal
        writeCoverageStats(ffh, fields, sample, famSample, famAnno)
        writeCoverageStats(gfh, fields, sample, geneSample, geneAnno)
        writeAlignerComparisonStats(afh, aFields, sample, alignerAgreement, geneSample['Total'])

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
    
        for k, v in alignerAgreement.iteritems():
            all_alignerAgreement[k] += v

    #Average
    n = len(sortedSamples)
    if n > 0:
        ffh.write("Average\t%.2f" %( famAll['Total']/n ))
        gfh.write("Average\t%.2f" %( geneAll['Total']/n ))
        for field in fields:
            if re.match("Anno", field):
                field = field.lstrip("Anno")
                ffh.write( "\t%.2f\t%.2f" %(famAllAnno[field]/n, famAllPcAnno[field]/n) )
                gfh.write( "\t%.2f\t%.2f" %(geneAllAnno[field]/n, geneAllPcAnno[field]/n) )
            else:
                ffh.write( "\t%.2f\t%.2f" %(famAll[field]/n, famAllPc[field]/n) )
                gfh.write( "\t%.2f\t%.2f" %(geneAll[field]/n, geneAllPc[field]/n) )
        ffh.write("\n")
        gfh.write("\n")

        #afh.write("Average")

    ffh.close()
    gfh.close()
    ifh.close()
    afh.close()

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
    for sample in sortedSamples:
        if sample not in sample2fam2gene2beds:
            sortedSamples.remove(sample)

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
        fam2gene2beds = sample2fam2gene2beds[sample] #genes mapped to Ref.
        famSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        geneSample = {'Total':0.0, 'Mapped':0.0, 'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
        
        famAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        geneAnno = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
        
        famSample['Total'] = len(fam2genes)
        famSample['Mapped'] = len(fam2gene2beds)
        
        geneSampleTotal = 0.0
        for fam, genes in fam2genes.iteritems():
            geneSampleTotal += len(genes)

        for fam, gene2beds in fam2gene2beds.iteritems():#each family with gene(s) mapped to Ref.
            assert fam in fam2genes
            #total = len(fam2genes[fam])
            mapped = len(gene2beds.keys())
            currAnnoCount = {'Insertion':0.0, 'Folded':0.0, 'Rearrangement':0.0}
            currCount = {'FullyMapped':0.0, 'Deletion': 0.0, 'Perfect':0.0, 'Insertion':0.0, 'Rearrangement':0.0, 'Folded':0.0}
            
            for gene, beds in gene2beds.iteritems():
                #if gene not in gene2fam or gene2fam[gene] not in fam2genes or gene in gene2status:
                if gene not in gene2fam or gene2fam[gene] not in fam2genes:
                    sys.stderr.write("Sample %s, gene %s is not in gene2fam or the geneList.\n" %(sample, gene))
                    #total -= 1
                    geneSampleTotal -= 1
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
            
            #geneSample['Total'] += total
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
        geneSample['Total'] = geneSampleTotal
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

#==========================================================================
#============================ CORE BASES ==================================
#==========================================================================
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

    return coreBases
#======================= END CORE BASES ==================================

#=============================================================
#====================== OPERON ===============================
#=============================================================
def writeOperon(f2, operon, genes, gene2beds, numPerfect, numRearrangement, numLost):
    f2.write("#%s\t%s\t%d\t%d\t%d\t%d\n" %(operon, ",".join(genes), len(genes), numPerfect, numRearrangement, numLost))
    for gene in genes:
        if gene not in gene2beds:
            f2.write("#%s\tLost\n" %gene)
            continue
        beds = gene2beds[gene]
        for bed in beds:
            f2.write("%s\t%d\t%d\t%s;%s\n" %(bed.chr, bed.start, bed.end, gene, operon))
    return

def checkOrderAndOrientation(genes, gene2beds):
    #Thu May  2 16:02:36 PDT 2013: 
    #!!! HACK: ALL WE CAN DO RIGHT NOW IS CHECK FOR ONLY THE ORDER OF THE GENES, 
    #          NOT THE ORIENTATION UNTIL halLiftover add in additional fields
    if len(genes) <= 1:
        return True
    
    reserved = True

    forward = True #the order of the genes has the same direction on Ref. + strand and on the original genome + strand
    beds1 = gene2beds[ genes[0] ] 
    beds2 = gene2beds[ genes[1] ]
    if beds2[0].start <= beds1[0].start:
        forward = False

    for i in xrange(0, len(genes) -1):
        beds1 = gene2beds[ genes[i] ]
        beds2 = gene2beds[ genes[i+1] ]
        if (forward and beds1[0].start > beds2[0].start) or \
           (not forward and beds2[0].start > beds1[0].start):
            return False
         
    return reserved

def checkOperons(operons, gene2beds, gene2len, g2status, outbasename):
    '''Check to see if: a/ all genes are reserved on C.Ref and 
                        b/ their orders and orientations are reserved on C.Ref
    '''
    geneReserved = 0 #number of operons that have all genes reserved on C.Ref
    ooReserved = 0 #number of operons that have the order and orientation of their genes reserved on C.Ref.
    f = open( "%s-operonStats.txt" %outbasename, "w")
    f.write("#OperonName\tgeneReserved\tooReserved\n") 
    f2 = open("%s-operons-broken.txt" %outbasename, "w")

    for operon in operons: 
        #print operon
        genes, strand = operons[operon] #list of genes belong to the operon and operon's strand
        #How many genes of the operon are reserved on C.Ref
        numPerfect = 0
        numRearrangement = 0
        numLost = 0
        gR = "No"
        ooR = "No"
        for gene in genes:
            assert gene in gene2len
            if gene not in gene2beds:
                numLost += 1
                continue
            status = getMapStatus(gene2beds[gene], gene2len[gene])
            #status = {'deletion':0, 'rearrangement':False, 'insertion':0, 'folded':False}
            if status['deletion'] == 0 and status['insertion'] == 0 and not status['rearrangement'] and not status['folded']:
                numPerfect += 1
            elif status['rearrangement'] or status['folded']:
                numRearrangement += 1

        if numPerfect == len(genes):
            geneReserved += 1
            gR = "Yes"
        if numLost == 0 and numRearrangement == 0:
            if checkOrderAndOrientation(genes, gene2beds):
                ooReserved += 1
                ooR = "Yes"
        f.write("%s\t%s\t%s\n" %(operon, gR, ooR))

        #if gR == "No" or ooR == "No":
        if ooR == "No":
            writeOperon(f2, operon, genes, gene2beds, numPerfect, numRearrangement, numLost)
    
    #Write summary stats
    f.write("\n#Total number of operons: %d\n" %len(operons))
    f.write("#Operons with all genes perfectly mapped: %d\t%.2f\n" %(geneReserved, getPc(geneReserved, len(operons)) ))
    f.write("#Operons with gene order and orientation reserved: %d\t%.2f\n" %(ooReserved, getPc(ooReserved, len(operons)) ))
    f.close()
    f2.close()

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

def readSampleBedFiles( samplepath ):
    if os.path.isdir(samplepath):
        files = [os.path.join(samplepath, file) for file in os.listdir(samplepath)]
    else:
        files = [samplepath]
    
    name2beds = {}
    for file in files:
        beds = readBedFile(file)
        for name, bedlist in beds.iteritems():
            if name in name2beds:
                name2beds[name].extend(bedlist)
            else:
                name2beds[name] = bedlist
    return name2beds

def readBedFiles(indir):
    sample2beds = {}
    for sample in os.listdir(indir):
        samplepath = os.path.join(indir, sample)
        if not os.path.isdir(samplepath):
            sample = sample.split('.')[0]
        sample2beds[sample] = readSampleBedFiles( samplepath )
    return sample2beds

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

def readGeneClusters2(file):
    gene2fam = {}
    fam2sample2genes = {} 

    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        if len(items) < 2:
            continue
        famid = items[0]
        sample2genes = {}
        for item in items[1:]:
            fields = item.split(';')
            sample = fields[0]
            genenames = fields[1].split(',')
            genes = []
            for g in genenames:
                #gi|386604534|ref|YP_006110834.1|
                geneitems = g.split('|')
                if len(geneitems) >= 4:
                    g = geneitems[3]
                genes.append(g)
            #for g in genes:
            #    gene2fam[g] = famid
            sample2genes[sample] = genes

        #Filter out family with duplicated genes:
        hasDups = False
        for s, genes in sample2genes.iteritems():
            if len(genes) > 1:
                hasDups = True
                break
        if not hasDups:
            fam2sample2genes[famid] = sample2genes
            for s, genes in sample2genes.iteritems():
                for g in genes:
                    gene2fam[g] = famid
    f.close()
    
    return gene2fam

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
    species = f.readline().strip().lstrip("#")
    for line in f:
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        assert len(items) >= 3
        name = items[0]
        genes = items[1].split(',')
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
                #sys.stderr.write("gene %s is not found in the geneClusters file. Sample %s\n" %(gene, sample))
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
                #sys.stderr.write("gene %s is not found in the geneClusters file. Sample %s\n" %(gene, sample))
                continue
            fam = gene2fam[gene]
            if fam not in fam2gene2beds:
                fam2gene2beds[fam] = {gene: beds}
            else:
                fam2gene2beds[fam][gene] = beds
        sam2fam2gene2beds[sample] = fam2gene2beds
    return sam2fam2gene2beds

def main():
    #bed directory: bed files, each containing genes of each sample mapped to Ref.
    #gene list directory: genelist files, each containing genes of each sample. Format: each line contains a <gene/protein name>\t<gene/protein sequence length>
    #geneClusters: file clusters genes into family. Format: <GeneFamilyId>\t<Comma,sep,list,of,genes,within,this,cluster>
    
    usage = "%prog <bed directory> <gene list directory> <output basename> <geneClusters> [<gene2sample psl directory>]"
    parser = OptionParser(usage = usage)
    parser.add_option('-r', '--ref', dest='ref', help='reference genome')
    parser.add_option('--refgenes', dest='refbed', help='Bed formatted file containing gene info of the reference genome.')
    parser.add_option('-s', '--geneStatus', dest='gstatus', help='Specify the list of annotated genes that were folded/rearranged')
    parser.add_option('-w', '--wiggle', dest='wiggle', help='Wiggle file showing the coverage of each position')
    parser.add_option('-o', '--operons', dest='operonFile', help='File containing operon set. Default=%default' )

    libplot.initOptions(parser)
    options, args = parser.parse_args()
    libplot.checkOptions(options, parser)

    if len(args) < 4:
        parser.error("Required at least 4 inputs\n")

    #check for genes that were misannotated (folded genes/ genes with indels)
    g2status = {}
    if options.gstatus:
        g2status = readGeneStatus(options.gstatus)
   
    #read input files
    sample2beds = readBedFiles(args[0]) #genes of each sample mapped to the Reference 
    sys.stderr.write("Done reading input bed files\n")
    sample2genes = readGeneLists(args[1]) #list of annotated genes of each sample
    sys.stderr.write("Done reading input gene lists\n")

    #number of coding bases share by all samples
    numsam = len( sample2beds.keys() )
    if options.wiggle:
        numsam2bases = getNumsam2codingBases(sample2beds, options.wiggle)
        for numsam in sorted( numsam2bases.keys() ):
            print "%d\t%d" %(numsam, numsam2bases[numsam])
    
    #cluster genes into gene families:
    #gene2fam = readGeneClusters2(args[3])  #ignore gene families with paralogous genes
    gene2fam = readGeneClusters(args[3])
    sam2fam2gene2beds = getSample2fam2gene2beds(sample2beds, gene2fam) 
    sam2fam2genes = getSample2fam2genes(sample2genes,gene2fam)

    if options.ref:
        assert options.refbed and  os.path.exists( options.refbed )
        ref_gene2beds = readSampleBedFiles( options.refbed )
        sys.stderr.write("Done reading reference annotated genes\n")
        coverage2(options.ref, ref_gene2beds, sam2fam2gene2beds, sam2fam2genes, gene2fam, g2status, args[2])
    else:
        coverage(sam2fam2gene2beds, sam2fam2genes, gene2fam, g2status, args[2])
    
    ##sharedGenes(sam2fam2gene2beds, sam2fam2genes, gene2fam, args[2], options) 
    
    if options.operonFile:
        sample, operons = readOperonFile(options.operonFile)
        checkOperons(operons, sample2beds[sample], sample2genes[sample], g2status, args[2])

if __name__ == "__main__":
    from referenceViz.src.geneStats import *
    main()






