#!/usr/bin/env python

import re, os, sys
from optparse import OptionParser
import numpy as np

#Sample	Total	NoDel	NoFolded	NoIns	NoP	NoRearrange	YesDelIA	YesDelO	YesFoldedIA	YesFoldedO	YesInsIA	YesInsO	YesNo	YesPA	YesPIA	YesPO	YesRearrangeIA	YesRearrangeO	YesYes	NoYes
def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        items.append(line)
    f.close()
    return items

def readFile(file):
    sample2cat2count = {}
    sample2total = {}
    f = open(file, 'r')
    fields = f.readline().strip().split('\t')

    for line in f:
        items = line.strip().split('\t')
        sample = items[0]
        sample2cat2count[sample] = {}
        for i in xrange(1, len(items)):
            category = fields[i]
            count = int(items[i])
            sample2cat2count[sample][category] = count
        total = sample2cat2count[sample]['Total']
        sample2total[sample] = total
        #for cat, count in sample2cat2count[sample].iteritems():
        #    sample2cat2count[sample][cat] = 100.0*count/total
    f.close()
    return sample2cat2count, sample2total

def readFiles(indir):
    sample2stats = {}
    stats = {} #sample1ToSample2ToCategory2Count
    sample2total = {}
    for ref in os.listdir(indir):
        refpath = os.path.join(indir, ref)
        if not os.path.isdir( refpath ):
            continue
        reffile = os.path.join(refpath, "%s-alignerComparisons.txt" %ref)
        sample2cat2count, currsample2total = readFile(reffile)
        for s, t in currsample2total.iteritems():
            if s not in sample2total:
                sample2total[s] = t
        stats[ref] = sample2cat2count
        for sample, cat2count in sample2cat2count.iteritems():
            if sample not in sample2stats:
                sample2stats[sample] = {}
            for cat, count in cat2count.iteritems():
                if cat not in sample2stats[sample]:
                    sample2stats[sample][cat] = [count]
                else:
                    sample2stats[sample][cat].append(count)
    return sample2stats, stats, sample2total

def calcAverage(sample2stats, outfile):
    f = open(outfile, 'w')
    fields = sorted(sample2stats.values()[0].keys())
    f.write("Sample\t%s\n" %"\t".join(fields))
    for sample, cat2counts in sample2stats.iteritems():
        f.write("%s" %sample)
        means = [ np.mean(cat2counts[field]) for field in fields ]
        #stds = [ np.mean(cat2counts[field]) for field in fields ]
        for m in means:
            f.write("\t%.2f" %m)
        f.write("\n")
    f.close()

def addLists(lists):
    sumList = []
    for i in xrange(0, len(lists[0])):
        count = [l[i] for l in lists]
        sumList.append( sum(count) )
    return sumList

def calcAverage2(sample2stats, outfile):
    aggFields = {1:'PerfectAligned', 2:'ImperfectAligned', 3:'Unaligned', 4:'UnalignedAligned', \
                 5: 'AlignedUnaligned', 6:'AlignedDisagreement', 7:'Total', 8:'Aligned1', 9:'Aligned2', \
                 10: 'PerfectAligned/Aligned1', 11:'Aligned/Aligned1'}
    
    #Sample	Total	NoDel	NoFolded	NoIns	NoP	NoRearrange	YesDelIA	YesDelO	YesFoldedIA	YesFoldedO	YesInsIA	YesInsO	YesNo	YesPA	YesPIA	YesPO	YesRearrangeIA	YesRearrangeO	YesYes	NoYes
    f = open(outfile, 'w')

    totalMeanCounts = {}
    for k in aggFields:
        totalMeanCounts[k] = 0.0
    
    f.write("#Sample")
    sortedFields = sorted(aggFields.keys())
    for i in sortedFields:
        f.write("\t%s" %aggFields[i])
    f.write("\n")

    for sample, cat2counts in sample2stats.iteritems():
        counts = {}
        for k in aggFields:
            counts[k] = 0.0

        #Total:
        counts[7] = cat2counts['Total']

        #PerfectAligned = YesPA
        counts[1] = cat2counts['YesPA']
        
        #ImperfectAligned = Yes*IA
        counts2 = [ cat2counts['YesPIA'], cat2counts['YesInsIA'], cat2counts['YesDelIA'], cat2counts['YesRearrangeIA'], cat2counts['YesFoldedIA'] ]
        counts[2] = addLists(counts2)
         
        #Unaligned = Total - (YesYes + YesNo + NoYes)
        yeslists = addLists( [ cat2counts['YesYes'], cat2counts['YesNo'], cat2counts['NoYes'] ] )
        counts[3] = [ counts[7][i] - yeslists[i] for i in xrange(0, len(counts[7])) ]

        #UnalignedAligned = No*
        counts4 = [ cat2counts['NoP'], cat2counts['NoIns'], cat2counts['NoDel'], cat2counts['NoRearrange'], cat2counts['NoFolded'] ]
        counts[4] = addLists(counts4)

        #AlignedUnaligned = YesNo
        counts[5] = cat2counts['YesNo']

        #AlignedDisagreement = Yes*O
        counts6 = [ cat2counts['YesPO'], cat2counts['YesInsO'], cat2counts['YesDelO'], cat2counts['YesRearrangeO'], cat2counts['YesFoldedO'] ]
        counts[6] = addLists(counts6)

        #Aligned1 = PerfectAligned + ImperfectAligned + Aligned*
        counts[8] = addLists( [counts[1], counts[2], counts[5], counts[6]] )

        #Aligned2 = PerfectAligned + ImperfectAligned + UnalignedAligned + AlignedDisagreement
        counts[9] = addLists( [counts[1], counts[2], counts[4], counts[6]] )

        counts[10] = [ float(counts[1][i])/counts[8][i] for i in xrange(0, len(counts[1])) ]
        counts[11] = [ float(counts[1][i] + counts[2][i])/counts[8][i] for i in xrange(0, len(counts[1])) ]

        meanCounts = {}
        for k, c in counts.iteritems():
            meanCount = np.mean(c) 
            meanCounts[k] = meanCount
            totalMeanCounts[k] += meanCount
         
        f.write("%s" %sample)
        for field in sortedFields:
            c = meanCounts[field]
            f.write("\t%.2f" %c)
        f.write("\n")
    
    f.write("Average")
    for field in sortedFields:
        totalCount = totalMeanCounts[field]
        avrCount = totalCount/len(sample2stats)
        f.write("\t%.2f" %avrCount)
    f.write("\n")

#=== get pairwise stat table ===
#Columns: <sample1>\t<sample2>\t<Total1>\t<Total2>\t<pairwise alignment shared genes>\t<cactus shared genes>
def getPc(count, total):
    if total == 0:
        return 0
    return 100.0*count/total

def getAlignedCounts(cat2count):
    bothAligned = cat2count['YesPA'] + cat2count['YesPIA'] + cat2count['YesInsIA'] + cat2count['YesDelIA'] + cat2count['YesRearrangeIA'] + cat2count['YesFoldedIA']
    alignedDisagreement = cat2count['YesPO'] + cat2count['YesInsO'] + cat2count['YesDelO'] + cat2count['YesRearrangeO'] + cat2count['YesFoldedO']
    aligned1 = bothAligned + alignedDisagreement + cat2count['YesNo'] 
    aligned2 = bothAligned
    #aligned2 = bothAligned + cat2count['NoP'] + cat2count['YesPO']
    return aligned1, aligned2 

def getPairwiseTab(stats, sample2total, file):
    f = open(file, 'w')
    f.write("#Sample1\tSample2\tTotal1\tTotal2\tPairwiseAligned\tProgressiveCactus\tProgressiveCactus/PairwiseAligned\n")
    for s1 in stats:
        total1 = sample2total[s1]
        for s2 in stats[s1]:
            total2 = sample2total[s2]
            cat2count = stats[s1][s2]
            align1, align2 = getAlignedCounts(cat2count)
            pc = getPc(align2, align1)
            f.write("%s\t%s\t%s\t%s\t%d\t%d\t%f\n" %(s1, s2, total1, total2, align1, align2, pc))
    f.close()

#============ MAIN =================
def main():
    usage = "%prog <indir> <outbasename> [options]"
    parser = OptionParser(usage = usage)
    #parser.add_option('--order', dest='order', help='the order to process the samples')
    options, args = parser.parse_args()
    #options.samples = []
    #if options.order:
    #    options.samples = readList(options.order)

    indir = args[0]
    sample2stats, stats, sample2total = readFiles(indir)
    outbase = args[1]
    #calcAverage(sample2stats, outfile)
    calcAverage2(sample2stats, outbase)

    getPairwiseTab(stats, sample2total, "%s-pairwise.tab"  %outbase)

if __name__ == '__main__':
    main()


