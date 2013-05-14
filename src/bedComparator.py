#!/usr/bin/env python

#nknguyen soe ucsc edu
#Mon May 13 11:32:29 PDT 2013
#
#Comparing two input alignments (each is represented by a bed-formatted file)
#Compute how much agreement/disagreement the two alignments have
#(The specific interest: gene alignments. For example: 
#comparing the lift-overs gene annotations of the query genome to the target 
#genome using aligner #1 versus lift-over using aligner #2)
#
#Input: 1/ Bed file 1 (alignment 1)
#       2/ Bed file 2 (alignment 2, with the same query & target with alignement 1)
#       3/ Optional: list of item/region names of interest (e.g list of gene names)
#
#Output: calculate counts of the following categories:
#       1/ Perfect Aligned Argreement: how many items (genes) where 
#                                      alignment 1 and 2 are identical
#       2/ Imperfect Aligned Agreement: alignment 1 and 2 have >= XX% overlap
#       3/ Unaligned Agreement: how many items of the query do not aligned to the
#                               target by both aligners
#       4/ Unaligned-Aligned Disagreement: items of query unaligned to target by
#                               aligner 1 but aligned to target by aligner 2
#       5/ Aligned-Unaligned Disagreement: items of query aligned to target by 
#                               aligner 1 but unaligned to target by aligner 2
#       6/ Aligned Disagreement: items of query aligned to target by both aligners
#                               but to different locations
#

import os, sys, re, time
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
        elif self.chromStart != other.chromStart:
            return cmp(self.chromStart, other.chromStart)
        else:
            return cmp(self.chromEnd, other.chromEnd)

############### ERROR CLASSES #################
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

############### ERROR CLASSES #################
class BedFormatError(Exception):
    pass

############### FUNCTIONS #############
#============ Read input files ================
def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        items.append( line.strip() )
    f.close()
    return list( set(items) )

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

def readMapfile(file):
    otherNames = {} #key = query name of genome1, val = list of query names of genome2 mapped to key
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        items = line.split('\t')
        assert len(items) >= 2
        names1 = items[0].split(',')
        names2 = items[1].split(',')
        for name in names1:
            if name not in otherNames:
                otherNames[name] = names2
            else:
                otherNames[name].extend(names2)
        for name in names2:
            if name not in otherNames:
                otherNames[name] = names1
            else:
                otherNames[name].extend(names1)
    f.close()
    return otherNames

#============ misc. functions ================
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

#============ comparing functions ================
def compareAlignments12( beds1, beds2, coverage ):
    #compare 2 alignments, beds have blocks info
    pass

def getNonOverlapRegions(beds):
    chr2regs = {}
    beds = sorted(beds)
    for bed in beds:
        reg = Region(bed.chrom, bed.chromStart, bed.chromEnd)
        if reg.chr not in chr2regs: #new chromosome
            chr2regs[reg.chr] = [reg]
        else:
            for r in chr2regs[reg.chr]:
                if r.start < reg.end and reg.start < r.end: #overlap
                    r.start = min( r.start, reg.start )
                    r.end = max( r.end, reg.end )
                    break
                else: #non-overlap
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
    chr2regs1 = getNonOverlapRegions(beds1) 
    chr2regs2 = getNonOverlapRegions(beds2) 
    total1 = 0.0
    total2 = 0.0
    overlap = 0.0
    for chr, regs1 in chr2regs1.iteritems():
        total1 += sum([reg.end - reg.start for reg in regs1])
        if chr not in chr2regs2:
            continue
        regs2 = chr2regs2[chr]
        total2 += sum([reg.end - reg.start for reg in regs2])
        overlap += calcOverlap(regs1, regs2) 
   
    pc1 = getPc(overlap, total1)
    pc2 = getPc(overlap, total2)
    if pc1 == 1.0 and pc2 == 1.0: 
        return 1 #perfect aligned
    elif pc1 >= coverage and pc2 >= coverage:
        return 2 #imperfect aligned
    else:
        return 6 #other

def printStats(categories, counts, outname):
    file = "%s.txt" %outname
    f = open(file, 'w')
    
    cats = sorted( counts.keys() )
    headers = [categories[i] for i in cats]
    f.write("#Stats\t%s\tTotal\tAligned1\tAligned2\n" %("\t".join(headers)))
    f.write("Count")
    for i in cats:
        f.write("\t%d" %counts[i])

    total = sum( counts.values() )
    aligned1 = total - counts[3] - counts[4]
    aligned2 = total - counts[3] - counts[5]
    f.write("\t%d\t%d\t%d\n" %(total, aligned1, aligned2))

    f.close()

def compareBeds(beds1, beds2, names, coverage, outname):
    '''query and target of beds1 and beds2 must be the same
    '''
    categories = {1:"PerfectAligned", 2:"ImperfectAligned", 3:"Unaligned", 4:"UnalignedAligned", 5:"AlignedUnaligned", 6:"AlignedDisagreement"}
    #Initialize counts (key = category, val = count)
    counts = {}
    for c in categories:
        counts[c] = 0

    for name in names:
        if name in beds1 and name not in beds2:
            counts[5] += 1
        elif name not in beds1 and name in beds2:
            counts[4] += 1
        elif name not in beds1 and name not in beds2:
            counts[3] += 1
        else:
            b1 = beds1[name]
            b2 = beds2[name]
            if b1[0].bed12:
                category = compareAlignments12(b1, b2, coverage)
            else:
                category = compareAlignments(b1, b2, coverage)
            counts[category] += 1
    
    printStats(categories, counts, outname)

#============ Case where query of beds1 != query of beds2, need to do the mapping ==========
def compareBeds2(beds1, beds2, names, coverage, outname, mapfile):
    '''target of beds1 and beds2 must be the same. query of beds1 and beds2 are mapped by mapfile
       In this case, report everything with respect to query of beds1.
    '''
    #map queries of beds1 and beds2:
    nameDict = readMapfile(mapfile)
    names2 = getOtherList(beds2.keys(), nameDict)#list of items in beds2 mapped to names of items in beds1
    #get union list of names:
    if not names:
        names = getUnionList( beds1.keys(), names2 ) 

    categories = {1:"PerfectAligned", 2:"ImperfectAligned", 3:"Unaligned", 4:"UnalignedAligned", 5:"AlignedUnaligned", 6:"AlignedDisagreement"}
    #Initialize counts (key = category, val = count)
    counts = {}
    for c in categories:
        counts[c] = 0

    for name in names:
        if name in beds1 and name not in names2:
            counts[5] += 1
        elif name not in beds1 and name in names2:
            counts[4] += 1
        elif name not in beds1 and name not in names2:
            counts[3] += 1
        else:
            b1 = beds1[name]
            currNames2 = nameDict[name]
            currCats = []
            for n2 in currNames2:
                b2 = beds2[n2]
                if b1[0].bed12:
                    category = compareAlignments12(b1, b2, coverage)
                else:
                    category = compareAlignments(b1, b2, coverage)
                currCats.append(category)
            #pick the best category 1 > 2 > 6
            counts[ min(currCats) ] += 1
    
    printStats(categories, counts, outname)

###################### MAIN ##################
def addOptions(parser):
    parser.add_option('-o', '--output', dest='outname', default='out', help='Output basename. Default=%default')
    parser.add_option('-n', '--names', dest='namelist', help='File with list of names of items of interest. If not specified, use the union list of items of the two input alignments.')
    parser.add_option('--coverage', dest='coverage', type='int', default=0.9, help='Minimum coverage for the "imperfect" aligned agreement category. Values: 0-1. Default=%default.')
    parser.add_option('--mapfile', dest='mapfile', help='If the queries of bed1 and bed2 are different, this file provides the mapping between queries of bed1 to queries of bed2 and vice versa. Format: each line represents 1 cluster (genes that mapped to each other): <gene1.1,gene1.2,...>\\t<gene2.1,gene2.2,...>. Default=%default')

def checkOptions(parser, options, args):
    if len(args) < 2:
        parser.error("Two input files are required. Only %d was given." %len(args))
    if not os.path.exists(args[0]):
        parser.error("Input file %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        parser.error("Input file %s does not exist.\n" %args[1])
    options.names = None
    if options.namelist:
        if not os.path.exists(options.namelist):
            parser.error("Namelist file %s does not exist.\n" %(options.namelist))
        else:
            options.names = readList(options.namelist)  

def main():
    usage = "%prog <bed1> <bed2>"
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, options, args)

    #Read input bed files:
    beds1 = readBedFile(args[0])
    beds2 = readBedFile(args[1])
    
    if not options.mapfile:
        if not options.names:
            options.names = getUnionList( beds1.keys(), beds2.keys() )
        compareBeds( beds1, beds2, options.names, options.coverage, options.outname )
    else:
        compareBeds2( beds1, beds2, options.names, options.coverage, options.outname, options.mapfile ) 

if __name__ == '__main__':
    main()

