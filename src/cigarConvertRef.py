#!/usr/bin/env python

'''
nknguyen soe ucsc edu
March 21 2012
Input: 1/  cigar file
       2/ list of reference names to convert to
Output: cigar file with coordinates relative to the reference in the list
'''

import os, sys, re
from optparse import OptionParser

class Ref():
    def __init__(self, line):
        line = line.strip()
        items = line.split('.')
        self.longname = line
        if len(items) == 6:
            self.name = '.'.join(items[0:2])
            self.chrlen = int(items[2])
            self.start = int(items[3])
            self.fraglen = int(items[4])

            if not isinstance(self.chrlen, int) or not isinstance(self.start, int) or not isinstance(self.fraglen, int):
                raise ValueError("Wrong reference header format. Should be: spc.chr.chrlen.start.fraglen.strand, where chrlen, start and fraglen are intergers. Got: %s\n" %line)
            self.strand = items[5]
            if self.strand != '1':
                raise ValueError("Got sequence header of negative strand: %s. Positive strand is required. Sorry.\n" %line)
        else:
            self.name = line
            self.chrlen = -1
            self.start = 0
            self.fraglen = -1
            self.strand = '.'

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

def getRef(frag, refs):
    start = frag.start
    end = frag.start + frag.fraglen

    for r in refs:
        if r.start <= start and end <= r.start + r.fraglen:
            return r
    #raise ValueError("Could not find sequence header that is a super-sequence of %s\n" %frag.longname)
    sys.stderr.write("%s\n" %frag.longname)
    return frag

def convertRef(name, start, end, refs):
    currRef = Ref(name)
    if currRef.name not in refs:
        raise ValueError("Ref %s was not found in the input ref file.\n" %currRef.name)
    #ref = refs[ currRef.name ]
    ref = getRef( currRef, refs[ currRef.name ] )
    newname = ref.longname
    newstart = start + currRef.start #Convert the coordinate to the original ref (genome coordinate)
    newend = end + currRef.start 
    if ref.chrlen > 0:#If the desired ref is only a subsequence of the original ref, convert the coordinate so that it is relative to this ref
        newstart -= ref.start
        newend -= ref.start

    if newstart < 0 or newend < 0:
        raise ValueError("Negative coordinates: old ref: %s, old start: %d, old end: %d; newref: %s, new start: %d, new end: %d\n" %(name, start, end, newname, newstart, newend))
    elif ref.fraglen > 0 and abs(newend - newstart) > ref.fraglen:
        raise ValueError("Out of range: old ref: %s, old start: %d, old end: %d; newref: %s, new start: %d, new end: %d\n" %(name, start, end, newname, newstart, newend))
    
    return newname, newstart, newend

def convertCigars(cigars, refs, outfile):
    f = open(outfile, 'w')
    for cigar in cigars:
        cigar.qname, cigar.qstart, cigar.qend = convertRef( cigar.qname, cigar.qstart, cigar.qend, refs )
        cigar.tname, cigar.tstart, cigar.tend = convertRef( cigar.tname, cigar.tstart, cigar.tend, refs )
        f.write("%s\n" %cigar.getStr())
    f.close()

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

def readReffile(file):
    f = open(file, 'r')
    refs = {} #key = 'shortName' (e.g: hg19.chr6), val = Ref object
    for line in f:
        ref = Ref(line)
        if ref.name in refs:
            refs[ref.name].append(ref)
            #raise ValueError("Repetitive reference %s will be ignored.\n" %ref.name)
        else:
            refs[ref.name] = [ref]
    f.close()
    return refs

def main():
    usage = ('%prog infile reffile outfile')
    parser = OptionParser( usage = usage )
    options, args = parser.parse_args()
    if len(args) < 3:
        raise ValueError("Required 3 files: inputfile, reffile, and outfile\n")
    infile = args[0]
    reffile = args[1]
    outfile = args[2]

    refs = readReffile(reffile)
    cigars = readCigarFile(infile)
    convertCigars(cigars, refs, outfile)

if __name__== "__main__":
    main()
