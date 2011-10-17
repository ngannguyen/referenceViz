#!/usr/bin/env python

import os, sys, re
import pysam

#Streaming bam file from stdin, print to stdout uniquelymapped stats

def decodeFlag(flag):
    flagDict = {"isPaired": 0, "properlyPaired":1, "unmapped":2, \
                "unmappedMate":3, "strand":4, "mateStrand":5, "firstRead":6,\
                "secondRead":7, "notPrimaryAlign":8, "failedQual":9,\
                "duplicate":10}

    for key in flagDict:
        if 2**flagDict[key] & flag == 0:
            flagDict[key] = False
        else:
            flagDict[key] = True
    return flagDict

#===========
f = pysam.Samfile( "-", "rb")

pair = 0
prevReadUniquelyMapped = False
prevFlagDict = {}
counts = {"uniquely_mapped":0,\
          "with_itself_and_mate_uniquely_mapped":0,\
          "with_itself_and_mate_uniquely_mapped_and_properly_paired":0,\
          "uniquely_mapped_mate_not_uniquely_mapped":0,\
          "with_itself_and_mate_not_uniquely_mapped":0}

#xacount = 0
prevRead = None

for r in f:
    pair += 1
    flagDict = decodeFlag( r.flag )
    uniquelyMapped = True
    if r.tags != None and 'XA' in [ tag for (tag, vals) in r.tags]:
        #xacount += 1
        uniquelyMapped = False
   
    #Make sure two reads are paired:
    if pair %2 == 0:
        if flagDict['isPaired'] and (r.qname != prevRead.qname or (r.is_read1 and prevRead.is_read1) or (r.is_read2 and prevRead.is_read2)):
            sys.stderr.write("Wrong pair! PrevRead: %s, Current read: %s\n" %(prevRead.qname, r.qname))
            sys.exit(1)
    #DEBUG:
        #if prevRead != None:
        #    if prevReadUniquelyMapped:
        #        sys.stderr.write("\nPrevRead: %s, uniquelyMapped\n" %(prevRead.qname))
        #    elif not prevFlagDict["unmapped"]:
        #        sys.stderr.write("\nPrevRead: %s, mapped, but not uniquely\n" %(prevRead.qname))
        #    else:
        #        sys.stderr.write("\nPrevRead: %s, Not mapped\n" %(prevRead.qname))
        #if uniquelyMapped and not flagDict["unmapped"]:
        #    sys.stderr.write("CurrentRead: %s, uniquelyMapped\n" %(r.qname))
        #elif not flagDict["unmapped"]:
        #    sys.stderr.write("CurrentRead: %s, mapped but not uniquely\n" %(r.qname))
        #else:
        #    sys.stderr.write("CurrentRead: %s, not mapped\n" %(r.qname))

    if uniquelyMapped:
        if not flagDict['unmapped']:#Current read is mapped uniquely
            counts["uniquely_mapped"] += 1
            
            if flagDict['isPaired'] and pair %2 == 0:
                if prevReadUniquelyMapped:
                    counts["with_itself_and_mate_uniquely_mapped"] += 2
                    #sys.stderr.write("with_itself_and_mate_uniquely_mapped += 2\n")
                    if flagDict['properlyPaired']:
                        counts["with_itself_and_mate_uniquely_mapped_and_properly_paired"] += 2
                elif not prevFlagDict['unmapped']:#prev read mapped but not uniquely
                    counts["uniquely_mapped_mate_not_uniquely_mapped"] += 2
                    #sys.stderr.write("uniquely_mapped_mate_not_uniquely_mapped += 2\n")
    else:
        if not flagDict["unmapped"]:#current read is mapped but Not uniquely
            if flagDict["isPaired"] and pair %2 == 0:
                if prevReadUniquelyMapped:
                    counts["uniquely_mapped_mate_not_uniquely_mapped"] += 2
                    #sys.stderr.write("uniquely_mapped_mate_not_uniquely_mapped += 2\n")
                elif not prevFlagDict['unmapped']:#prev read is mapped but not uniquely:
                    counts["with_itself_and_mate_not_uniquely_mapped"] +=2
                    #sys.stderr.write("with_itself_and_mate_not_uniquely_mapped += 2\n")

    #prevReadUniquelyMapped = uniquelyMapped and (not flagDict["unmapped"])
    if (not flagDict["unmapped"]) and uniquelyMapped:
        prevReadUniquelyMapped = True
    else:
        prevReadUniquelyMapped = False
    prevFlagDict = flagDict.copy()
    prevRead = r

#Print out stats:
#for k in sorted( counts.keys() ):
#    sys.stdout.write("%d\t%s\n" %(counts[k], k))
fields = ["uniquely_mapped", "with_itself_and_mate_uniquely_mapped", "with_itself_and_mate_uniquely_mapped_and_properly_paired", "uniquely_mapped_mate_not_uniquely_mapped", "with_itself_and_mate_not_uniquely_mapped" ]
for k in fields:
    sys.stdout.write("%d\t%s\n" %(counts[k], k))

#sys.stderr.write("Num XA reads: %d\n" %(xacount))

