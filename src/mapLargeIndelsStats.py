#!/usr/bin/env python

"""
Extract insertions (in each sample with respect to hg19) with size >= y bases
Get the location of the insertions on hg19, as well as the sequences being inserted.
Mapped the inserted sequences to hg19.
Record where those inserted sequences mapped on hg19 to find out the source and sink of
the duplication/rearragement events. 
Input: pathStats_hg19.xml
Output:
"""

import os, sys, re, copy
import xml.etree.ElementTree as ET
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import system
from sonLib.bioio import logger
from sonLib.bioio import setLogLevel
from sonLib.bioio import getTempDirectory

class Insertion():
    def __init__(self, line, sampleName, referenceName):
        items = line.strip().split()
        if len(items) < 16:
            sys.stderr.write("Cactus indel record (pathStats_*.xml) does not have enough fields\n")
            sys.exit(1)
        self.name = items[1]
        self.start = int( items[3] )
        self.length = int( items[5] )
        self.strand = items[7]
        self.ref = items[9]
        self.refstart = int(items[11])
        self.reflength = int(items[13])
        self.refstrand = items[15]
        refitems = self.ref.split('.')
        if len(refitems) < 2:
            sys.stderr.write("Reference sequence must provide chromosome info: i.e: species.chr\n")
            sys.stderr.write("%s\n" %self.ref)
            sys.exit(1)
        self.refspc = refitems[0]
        self.refchrom = refitems[1]
        #Convert coordinate if necessary:
        if self.refstrand == 0:
            sys.stder.write("Cactus indels are on negative strand\n")
            sys.exit(1)
        if len(refitems) == 6:
            refstrand = refitems[5]
            reflen = int(refitems[2])
            offset = int(refitems[3])
            fraglen = int(refitems[4])
            if refstrand == "1": #Positive strand
                self.refstart = self.refstart + offset
            else:
                self.refstart = self.refstart + reflen - (offset + fraglen)
        self.sampleName = sampleName
        self.referenceName = referenceName

    def setSeq(self, seq):
        self.seq = seq

    def __cmp__(self, other):
        chr = int( self.refchrom.lstrip('chr') )
        otherchr = int(other.refchrom.lstrip('chr'))
        if chr < otherchr:
            return -1
        elif chr > otherchr:
            return 1
        else:
            if self.refstart < other.refstart:
                return -1
            elif self.refstart == other.refstart:
                return 0
            else:
                return 1

class Cigar():
    def __init__(self, line):
        items = line.strip().split()
        if len(items) < 11:
            sys.stderr.write("Wrong cigar format, must have at least 11 fields: %s\n" %(line) )
            sys.exit(1)
        self.query = items[1]
        self.qstart = int(items[2])
        self.qend = int(items[3])
        self.qstrand = items[4]
        qitems = self.query.split('-')
        #id = "%s_%d_%d_%s_%d_%d_%s" %(ins.name, ins.start, ins.length, ins.strand, ins.refstart, ins.reflength, ins.refstrand)
        if len(qitems) < 7:
            sys.stderr.write("Wrong query name format, must have at least 7 fields: %s\n" %self.query)
            sys.exit(1)
        qlen = int(qitems[2])

        if self.qstrand == '-':
            qstart = qlen - self.qstart
            qend = qlen - self.qend
            self.qstart = qstart
            self.qend = qend
            #sys.stderr.write("Cigar, converted - starnd coord of %s to %d, %d\n" %(line, self.qstart, self.qend))
        self.target = items[5]
        self.tstart = int(items[6])
        self.tend = int(items[7])
        self.tstrand = items[8]
        if self.tstrand == '-':
            sys.stderr.write("Warning, tstrand is negative!! %s" %(line))
        self.score = int(items[9])
        self.cigarstr = " ".join(items[10:])

    def desc(self):
        return "%s %d %d %s %s %d %d %s %d %s" %(self.query, self.qstart, self.qend, self.qstrand, self.target, self.tstart, self.tend, self.tstrand, self.score, self.cigarstr)

    def __cmp__(self, other):#Compare in term of query coordinate
        if self.query < other.query:
            return -1
        elif self.query > other.query:
            return 1
        else:
            if self.qstart < other.qstart:
                return -1
            elif self.qstart > other.qstart:
                return 1
            else:
                if self.qend < other.qend:
                    return -1
                elif self.qend > other.qend:
                    return 1
                else:
                    return 0

def readfile(file, minLen, filteredSamples):
    xmltree = ET.parse(file) 
    root = xmltree.getroot()
    insertions = {} #val = sample, key = list of Insertions
    for sample in root.findall('statsForSample'):
        name = sample.attrib['sampleName']
        if name == 'reference' or name == 'aggregate' or name == 'ROOT' or name == '' or name in filteredSamples:
            continue
        #samples.append(name)
        insertions[name] = []
        ref = sample.attrib['referenceName']
        indels = sample.text.strip().split('\n')
        for site in indels:
            if site == '':
                continue
            indel = Insertion(site, name, ref)
            if indel.length >= minLen:
                insertions[name].append(indel)
    return insertions

def readFa(file):
    f = open(file, 'r')
    seqs = {} #key = header, val = sequence
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '':
                if header in seqs:
                    sys.stderr.write("Sequences with same headers? %s\n" %header)
                seqs[header] = seq
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if header != '' and seq != '':
        seqs[header] = seq
    return seqs

def readCigarFile( file ):
    cigars = []
    f = open(file, 'r')
    for line in f:
        cigars.append(Cigar(line))
    return cigars

def getRevComp(seq):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c'}
    revseq = ''
    for i in xrange( len(seq)-1, -1, -1 ):
        revseq += d[ seq[i] ]
    return revseq

def getSeq(seq, insertion):
    start = insertion.start
    if insertion.strand == '-':
        start = len(seq) - (insertion.length + start)
        seq = getRevComp(seq)
    end = start + insertion.length
    if start < 0 or len(seq) < end:
        sys.stderr.write("Indel with start %d, length %d, strand %s is out of range. Sequence %s has length %d\n" %(insertion.start, insertion.length, insertion.strand, insertion.name, len(seq)))
        sys.exit(1)
    return seq[start: end]

def printInsertSeqs(insertions):
    f = open('largeIndelSeqs.txt', 'w')
    for sample in insertions:
        ins = insertions[sample]
        for insertion in ins:
            f.write(">%s\n" %insertion.name)
            f.write("%s\n" %insertion.seq)

    f.close()

def getInsertedSeqs(insertions, seqdir):
#def getInsertedSeqs(insertions, insSeqFile):
    #seqs = readFa(insSeqFile)
    for sample in insertions:
        seqfile = os.path.join(seqdir, sample)
        if not os.path.exists(seqfile):
            sys.stderr.write("Sequence file %s is required, but does not exists\n" %(seqfile))
            sys.exit(1)
        seqs = readFa(seqfile) 
        ins = insertions[sample]
        for insertion in ins:
            seq = seqs[insertion.name]
            insertion.setSeq( getSeq(seq, insertion) )
            #insertion.setSeq( seq )
    return

def basesMappedDEBUG(cigars):
    bases = 0
    multiMapped = []
    currstart = 0
    if len(cigars) == 0:
        return bases, multiMapped
    cigars.sort()
    for cig in cigars:
        start = cig.qstart
        end = cig.qend
        print cig.desc()
        print "Start %d, end: %d, CurrStart: %d; bases: %d, multiMapped: %d" %(start, end, currstart, bases, len(multiMapped))
        if start < currstart: #Multimapping
            for pos in range(start, currstart):
                if pos not in multiMapped:
                    multiMapped.append(pos)
            start = currstart
        print "Start: %d, End: %d, Size: %d\n" %(start, end, end - start)
        if end > start:
            bases += end - start
        if currstart < end:
            currstart = end
    return (bases, len(multiMapped))

def basesMapped(cigars):
    bases = 0
    multiMapped = []
    currstart = 0
    if len(cigars) == 0:
        return bases, multiMapped
    cigars.sort()
    for cig in cigars:
        start = cig.qstart
        end = cig.qend
        if start < currstart: #Multimapping
            for pos in range(start, min([end, currstart])):
                if pos not in multiMapped:
                    multiMapped.append(pos)
            start = currstart
        if end > start:
            bases += end - start
        if currstart < end:
            currstart = end
    return (bases, len(multiMapped))

def getMhcMap(cigars, start, end):
    mhcCigs = []
    for cig in cigars:
        if cig.target == 'chr6' and (cig.tend >= start and cig.tstart <= end):
            mhcCigs.append(cig)
    return mhcCigs

def getBestMap(cigars, start, end):
    if len(cigars) == 0:
        return None
    bestcig = cigars[0]
    for cig in cigars:
        if cig.score > bestcig.score:
            bestcig = cig
        elif cig.score == bestcig.score and (cig.target == 'chr6' and cig.tend >= start and cig.tstart <=end):
            bestcig = cig
    return bestcig

def repeatPercentage( seq ):
    if len(seq) == 0:
        return 0.0
    numlower = 0.0
    for c in seq:
        if c == c.lower():
            numlower += 1
    return numlower*100/len(seq)

def multiCopyNumber(file, inslen, cutoff):
    #sys.stderr.write("GET COPY NUMBER, file %s!!\n" %file)
    cigars = readCigarFile(file)
    #sys.stderr.write("DONE\n")
    if len(cigars) == 0:
        return bases, multiMapped
    cigars.sort()

    multiMapped = []
    currstart = 0
    for cig in cigars:
        start = cig.qstart
        end = cig.qend
        #if file == "out3/mcf/selfAlign/mcf.chr6_mcf_hap5.4833398.2305647.406804.1-367196-4440-1-31292376-5016-1":
        #    sys.stderr.write("Currstart: %d\tStart %d\tend %d. LenMultiMapped: %d\n" %(currstart, start, end, len(multiMapped)))

        if start < currstart: #Multimapping
            for pos in range(start, min([end, currstart]) ):
                if pos not in multiMapped:
                    multiMapped.append(pos)
            start = currstart
        if currstart < end:
            currstart = end
            
    if len(multiMapped)*100.0/inslen >= cutoff:
        sys.stderr.write("MultipleMapping: %s\t%d, %d\n" %(file, len(multiMapped), inslen))
        return True
    return False
        
def stats( insertions, options ): 
    stats = {'Mapped':[0,0], 'Mapped To MHC':[0,0], 'Mapped Outside of MHC':[0,0], 'Tandem duplications':[0,0], 'Repetitive':[0,0], 'Copy Number Change':[0,0], 'Multimapping bases':[0,0]}
    label2id = {'Unmapped':[], 'Mapped':[], 'Mapped To MHC':[], 'Mapped Outside of MHC':[], 'Tandem duplications':[], 'Repetitive':[], 'Copy Number Change':[]}
    total = 0
    totallen = 0
    multiM = 0
    #mappedLen = 0

    statsFile = os.path.join(options.outdir, "stats-%d.txt" %options.cutoff)
    for sample in insertions:
        sampledir = os.path.join(options.outdir, sample)
        selfdir = os.path.join(sampledir, "selfAlign")

        insList = insertions[sample]
        for ins in insList:
            total += 1
            totallen += ins.length
            id = "%s-%d-%d-%s-%d-%d-%s" %(ins.name, ins.start, ins.length, ins.strand, ins.refstart, ins.reflength, ins.refstrand)
            sys.stderr.write("%s\n" %id)
            mapfile = os.path.join(sampledir, id)
            selfmapfile = os.path.join(selfdir, id)
            #multipleCopies = multiCopyNumber(selfmapfile, ins.length, self.options.cutoff)
            multipleCopies = multiCopyNumber(selfmapfile, ins.length, 80)
            if multipleCopies:
                stats['Copy Number Change'][0] += 1
                stats['Copy Number Change'][1] += ins.length
                label2id['Copy Number Change'].append(id)

            cigars = readCigarFile(mapfile)

            #Get percentage of total bases with lowercase
            rp = repeatPercentage( ins.seq )
            
            #Number of bases that got aligned, and number of bases that got aligned to multiple locations
            (mapped, multiMapped) = basesMapped(cigars)
            percentMapped = mapped*100.0/ins.length
            if percentMapped > 100:
                sys.stderr.write("Over counting number of mapped bases: %d, %d, %d%%\n" %(mapped, ins.length, percentMapped))
                basesMappedDEBUG(cigars)
                sys.exit(1)

            #Find cigar records that overlap with the MHC
            mhccigars = getMhcMap(cigars, options.start, options.end)
            mhcMapped = False
            #Those that mapped to the MHC, get cigar with best score:
            #bestCigar = getBestMap(mhccigars)
            #if not bestCigar: #None of the alignment mapped to the MHC
            #    bestCigar = getBestMap(cigars)
            #else:
            #    mhcMapped = True

            bestCigar = getBestMap(cigars, options.start, options.end)
            if bestCigar in mhccigars:
                mhcMapped = True

            if rp > 50:
                stats['Repetitive'][0] +=1
                stats['Repetitive'][1] +=ins.length
                label2id['Repetitive'].append(id)
            
            if not bestCigar or percentMapped < options.cutoff:
                label2id['Unmapped'].append(id)
            else:
                #multiM += multiMapped
                if multiMapped*100.0/ins.length >= options.cutoff: #More than 50% of the bases multimapped
                    stats['Multimapping bases'][0] += 1
                    stats['Multimapping bases'][1] += ins.length
                #bestCigar != None and percentMapped >= options.cutoff:
                stats['Mapped'][0] +=1
                stats['Mapped'][1] += ins.length
                #stats['Mapped'][1] += mapped
                label2id['Mapped'].append(id)
                #if rp > 50:
                #    stats['Repetitive'] +=1
                #    label2id['Repetitive'].append(id)
                if mhcMapped:
                    stats['Mapped To MHC'][0] +=1
                    stats['Mapped To MHC'][1] += ins.length
                    #stats['Mapped To MHC'][1] += mapped
                    label2id['Mapped To MHC'].append(id)
                else:
                    stats['Mapped Outside of MHC'][0] +=1
                    stats['Mapped Outside of MHC'][1] += ins.length
                    #stats['Mapped Outside of MHC'][1] += mapped
                    label2id['Mapped Outside of MHC'].append(id)
                if multipleCopies and ((ins.refstart >= bestCigar.tend and ins.refstart - bestCigar.tend <= options.dist) or ( bestCigar.tstart >= ins.refstart and bestCigar.tstart - ins.refstart <= options.dist)):
                    stats['Tandem duplications'][0] +=1
                    stats['Tandem duplications'][1] += ins.length
                    label2id['Tandem duplications'].append(id)

    f = open(statsFile, 'w')
    f.write("Total: %d, %d\n" %(total, totallen))
    for k in stats:
        f.write( "%s\t%d\t%d\t%.3f%%\n" %(k, stats[k][0], stats[k][1], stats[k][1]*100.0/totallen) )
        #if k in ['Tandem duplications', 'Repetitive', 'Copy Number Change']:
        #    f.write( "%s\t%d\t%.3f%%\n" %(k, stats[k], stats[k]*100.0/total) )
        #else:
        #    f.write( "%s\t%d\t%d\t%.3f%%\n" %(k, stats[k][0], stats[k][1], stats[k][1]*100.0/totallen) )
    #f.write("%s\t%d\t%.3f%%\n" %("Multimapping bases", multiM, multiM*100.0/totallen))
    f.write("\n")
    for k in label2id:
        f.write("%s\t%s\n" %(k, ','.join(label2id[k]) ))

    f.close()

def initOptions( parser ):
    parser.add_option('-f', '--filteredSamples', dest='filteredSamples', help='Comma separated list of samples that were filtered out (not to include in the plot)')
    parser.add_option('-s', '--minsize', dest='minsize', type='int', default = 1000, help='Minimum size of indels to be included in the analysis. Default = 1000')
    parser.add_option('-d', '--dist', dest='dist', type='int', default = 1000, help='Maximum distance to be called tandem duplication. Default = 1000')
    parser.add_option('-r', '--reference', dest='ref', default='/hive/users/nknguyen/mhc/data/assemblies/chroms/', help='reference to map the indels to')
    parser.add_option('-o', '--outdir', dest='outdir', default='.', help='Output directory. Default = current directory')
    parser.add_option('-i', '--identity', dest='identity', type='float', default=90, help='Percentage identity used as cutoff in alignments')
    parser.add_option('-c', '--cutoff', dest='cutoff', type='float', default=50, help='Percentage mapped used as cutoff in alignments')
    parser.add_option('--start', dest='start', type='int', default=28477754, help='Start of the MHC region. Default = 28477754')
    parser.add_option('--end', dest='end', type='int', default=33448354, help='Start of the MHC region. Default = 33448354')
    parser.add_option('--sequences', dest='seqs', default='/hive/users/benedict/referenceScripts/dataDir/mhcHumanVariantsNsRemovedAndFiltered/', help='Directory to the sample fasta sequences.')

def checkOptions( args, options, parser ):
    if len(args) < 2:
        parser.error('Please provide input pathStats.xml file and sequencesDir\n')
    if not os.path.exists(args[0]):
        parser.error('File %s does not exist\n' %args[0])
    if not os.path.exists(args[1]):
        parser.error('File %s does not exist\n' %args[1])
    if options.filteredSamples:
        options.filteredSamples = options.filteredSamples.split(',')
    else:
        options.filteredSamples = []
    if not os.path.exists(options.ref):
        parser.error("Directory %s does not exists\n" %options.ref)
    if os.path.isdir(options.ref): 
        options.refs = [ os.path.join(options.ref, d) for d in os.listdir(options.ref) ]
    else:
        options.refs = [ options.ref ]
    if not os.path.exists(options.outdir):
        system("mkdir -p %s" %options.outdir)

def main():
    usage = ("Usage: %prog [options] pathStats_hg19.xml sequencesDir")
    parser = OptionParser( usage = usage )
    initOptions(parser)
    options, args = parser.parse_args()
    checkOptions( args, options, parser )

    insertions = readfile(args[0], options.minsize, options.filteredSamples)
    #HLA-DRB region: chr6:32,414,165-32,605,002
    start = 32414165
    end = 32605002
    numDRBins = 0
    for s in insertions:
        for ins in insertions[s]:
            if start <= ins.refstart and ins.refstart <= end:
                numDRBins += 1
    sys.stderr.write("Number of insertions in HLA-DRB region: %d\n" %numDRBins)
    #return

    getInsertedSeqs(insertions, args[1])
    sys.stderr.write("Done getting insertion sequences\n")
    #printInsertSeqs(insertions)

    stats(insertions, options)

if __name__ == "__main__":
    main()






