#!/usr/bin/env python

import os, sys, re

class Cigar():
    def __init__(self, line):
        line = line.strip()
        items = line.split()
        if len(items) < 12:
            raise ValueError("Wrong Cigar format. At least 12 fields are expected. Only saw %d. Line: %s\n" %(len(items), line) )
        #query
        self.qname = items[1]
        self.qstart = int(items[2])
        self.qend = int(items[3])
        self.qstrand = items[4]

        if not isinstance(self.qstart, int) or not isinstance(self.qend, int):
            raise ValueError("Wrong Cigar format. Required start and end (field 3 and 4) to be integer. Got: %s, %s\n" %(items[2], items[3]))
        #target
        self.tname = items[5]
        self.tstart = int(items[6])
        self.tend = int(items[7])
        self.tstrand = items[8]
        self.score = items[9]
        self.cigarstr = " ".join(items[10:])

    def getStr(self):
        return "cigar: %s %d %d %s %s %d %d %s %s %s" %(self.qname, self.qstart, self.qend, self.qstrand, self.tname, self.tstart, self.tend, self.tstrand, self.score, self.cigarstr)

class Seq():
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        
        items = header.split(".")
        assert len(items) == 6
        subitems = items[1].split("-")
        assert len(subitems) >= 2
        self.chr = subitems[0]
        self.start = int( subitems[1] )
        length = int( items[4] )
        self.end = self.start + length

def readCigarFile(file):
    cigars = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        cigar = Cigar(line)
        cigars.append(cigar)
    f.close()
    return cigars

def readFastaFile(file):
    seqs = []
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == '>':
            if header != '' and seq != '':
                seqs.append( Seq(header, seq) )
            header = line.lstrip('>')
            seq = ''
        else:
            seq += line
    if seq != '' and header != '':
        assert header not in seqs
        seqs.append( Seq(header, seq) )
    f.close()
    return seqs

def isOverlap(seq1, seq2):
    if seq1.chr == seq2.chr and seq1.start < seq2.end and seq2.start < seq1.end:
        return True
    return False

####################################################
#Input directory contains cigar files and fasta files
#Check:
# For each set of (cigar, fasta):
#   1/the fasta contains all the target names in the cigar and
#   the cigar only contains target names in the fasta
#
#   2/no overlapping sequences in the fasta
#
#No sequence belong to two sets
#
indir = sys.argv[1]
files = os.listdir(indir)
basename = "-".join( files[0].split('.')[0].split('-')[:-1] )
set2files = {}
for file in files:
    items = file.split('.')
    setid = items[0].split('-')[-1]
    if setid not in set2files:
        set2files[setid] = [file]
    else:
        set2files[setid].append(file)

#Read input files
set2seqs = {}
set2cigs = {}
for id in set2files:
    cigfile = os.path.join(indir, "%s-%s.cig" %(basename, id))
    set2cigs[id] = readCigarFile(cigfile)

    fafile = os.path.join(indir, "%s-%s.fa" %(basename, id))
    set2seqs[id] = readFastaFile(fafile) 

error = 0
#Check the validity of each set:
for id, seqs in set2seqs.iteritems():
    cigs = set2cigs[id]
    tnames = []
    for cig in cigs:
        if cig.tname not in tnames:
            tnames.append(cig.tname)

    for s in seqs:
        if s.header not in tnames: #fafile contains contigs with not constraints
            error += 1
            print "set %s, sequence header %s is not in cigfile. CigHeaders: %s\n" %(id, s.header, ",,".join(tnames))

    headers = [s.header for s in seqs]
    for tname in tnames:
        if tname not in headers:
            error += 1
            print "set %s, cigar tname %s is not in contig fasta file. fasta headers: %s\n" %(id, tname, ",,".join(header))

    #check for overlapping contigs in the fasta file:
    for i, s1 in enumerate(seqs):
        for j, s2 in enumerate(seqs):
            if i == j:
                continue
            if s1.header == s2.header:
                error += 1
                print "Redundant contigs %s in fafile of set %s" %(s1.header, id)
            if isOverlap(s1, s2):
                error += 1
                print "Overlap contigs %s, %s in fafile of set %s" %(s1.header, s2.header, id)

#Make sure that no sequence belongs to more than 1 set
header2sets = {}
for id, seqs in set2seqs.iteritems():
    for s in seqs:
        if s.header not in header2sets:
            header2sets[s.header] = [id]
        else:
            error += 1
            print "Sequence %s belong to multiple sets %s, %s" %(s.header, id, ", ".join(header2sets[s.header]))
            header2sets[s.header].append(id)

if error == 0:
    print "Files are OK"



