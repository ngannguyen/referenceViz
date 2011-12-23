#!/usr/bin/env python

'''
10/19/2011: nknguyen soe ucsc edu
Take a union of the dbsnps snps & 1k snps:
'''

import os, sys, re, copy

class Snp():
    def __init__(self, line):
        self.desc = line
        items = line.split('\t')
        if len(items) < 25:
            sys.stderr.write("Each snp record must have at least 25 fields\n")
            sys.stderr.write("%s\n" %line)
        self.chrom = items[0]
        self.chromStart = int(items[1])
        self.chromEnd = int(items[2])
        self.name = items[3]
        self.score = int(items[4])
        self.strand = items[5]
        self.refNCBI = items[6].upper()
        self.refUCSC = items[7].upper()
        self.observed = items[8].upper().split('/')
        if self.strand == '-':
            self.observed = [ reverse(a) for a in self.observed] 
        #Check if indel is snp:
        self.isSnp = False
        for o in self.observed:
            if len(o) == self.chromEnd - self.chromStart and o!= '-' and o != self.refUCSC:
                self.isSnp = True
                break
        
        self.molType = items[9]
        self.type = items[10]
        #self.valid = items[11]
        self.func = items[14]
        self.locType = items[15]
        self.exceptions = items[17]
        self.alleleFreqCount = int(items[20])
        #self.alleles = items[21].rstrip(',')
        #self.alleleFreqs = items[23].rstrip(',')
        self.alleles = items[21].upper().rstrip(',').split(',')
        #for allele in self.alleles:
        #    if len(allele) < self.chromEnd - self.chromStart:
        #        sys.stderr.write("dbsnp line '%s', allele length is different from chromEnd - chromStart\n" %(line))
        #        sys.exit(1)
        if self.strand == '-':
            d = {'a':'t', 't':'a', 'c':'g', 'g':'c'}
            newalleles = []
            for a in self.alleles:
                if a in d:
                    newalleles.append( d[a] )
                else:
                    newalleles.append( a )
            self.alleles = newalleles
        #self.alleleFreqs = [ float(f) for f in items[23].rstrip(',').split(',') ]

    def __cmp__(self, other):
        chr = int( self.chrom.lstrip('chr') )
        otherchr = int( other.chrom.lstrip('chr') )

        if chr < otherchr:
            return -1
        elif chr > otherchr:
            return 1
        else:
            if self.chromStart < other.chromStart:
                return -1
            elif self.chromStart == other.chromStart:
                return 0
            else:
                return 1

class Snp2():
    def __init__(self, line):
        items = line.strip().split('\t')
        self.chrom = items[0]
        self.start = int(items[1]) -1
        self.ref = items[2].upper()
        self.alt = items[3].upper()
        self.filter = items[5]
        self.type = ''
        if len(self.ref) == 1 and len(self.alt) == 1:
            self.type = 'snp'
        elif self.alt != '<DEL>' and len(self.ref) > len(self.alt):
            self.type = 'deletion'
            self.start +=1
            self.ref = self.ref[1:]
            self.alt = '-'
        elif len(self.ref) < len(self.alt):
            self.type = 'insertion'
            self.ref = '-'
            self.start += 1
            self.alt = self.alt[1:]

def reverse(s):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    if s in d:
        return d[s]
    else:
        return s

def splitDbSnp(snp):
    #sys.stderr.write("%s" %snp.desc)
    snps = []
    #if len(snp.observed) == 0 or len(snp.refNCBI) != len(snp.observed[0]) or len(snp.refUCSC) != len(snp.observed[0]):
    if len(snp.observed) == 0:
        return snps
    snpLen = len(snp.observed[0])
    #sys.stderr.write("snplen %d\n" %snpLen)
    
    for allele in snp.observed:
        if len(allele) != snpLen:
            return snps
    for i in xrange(snpLen):
        s = copy.copy( snp )
        if len(snp.refNCBI) > i:
            s.refNCBI = snp.refNCBI[i]
        else:
            s.refNCBI = ''
        if len(snp.refUCSC) > i:
            s.refUCSC = snp.refUCSC[i]
        else:
            s.refUCSC = ''
        s.chromStart = snp.chromStart + i
        s.observed = [ a[i] for a in snp.observed ]
        s.isSnp = True
        s.type = 'single'
        snps.append(s)
    return snps

def readSnps(file, dbsnps):
    f = open(file, 'r')
    i = 0
    j = 0

    for line in f:
        snp = Snp2(line)
        if snp.filter != 'PASS':
            continue

        i = j
        if i >= len(dbsnps):
            sys.stderr.write("Index out of range: %d, %d\n" %(i, len(dbsnps)))
            sys.exit(1)
        dbsnp = dbsnps[i]
        flag = False
        while dbsnp.chromStart <= snp.start:
            if dbsnp.chromStart == snp.start:
                if snp.type == 'snp' and dbsnp.isSnp:
                    if snp.ref in dbsnp.observed and snp.alt in dbsnp.observed:
                        flag = True
                        break
                elif snp.type == 'deletion' and (dbsnp.type == 'deletion' or dbsnp.type == 'in-del'):
                    if len(snp.ref) - len(snp.alt) == dbsnp.chromEnd - dbsnp.chromStart:
                        flag = True
                        break
                elif snp.type == 'insertion' and (dbsnp.type == 'insertion' or dbsnp.type == 'in-del'):
                    inslen = len(snp.alt) - len(snp.ref)
                    for a in dbsnp.observed:
                        if a != '-' and len(a) == inslen:
                            flag = True
                            break

            if dbsnp.chromStart < snp.start:
                j += 1
            i += 1
            if i == len(dbsnps):
                break
            dbsnp = dbsnps[i]
        
        if not flag:
            end = 0
            type = ''
            observed = "/".join([snp.ref, snp.alt])
            if snp.type == 'snp':
                end = snp.start + 1
                type = 'single'
            elif snp.type == 'insertion':
                end = snp.start
                type = 'insertion'
            elif snp.type == 'deletion':
                end = snp.start + len(snp.ref)
                type = 'deletion'

            sys.stdout.write('chr%s\t%d\t%d\t\t0\t\t%s\t%s\t%s\t\t%s\t\t0\t0\t\t\t\t\t\t\t\t\t\t\t\t\n' %(snp.chrom, snp.start, end, snp.ref, snp.ref, observed,type))

    return

def readDbSnps(file):
    f = open(file, 'r')

    snps = []
    nummnp = 0

    for line in f:
        if re.search('chromStart', line):
            continue
        snp = Snp(line)
        #if snp.type == 'single' or snp.type == 'mnp':
        if snp.type == 'mnp':
            subsnps = splitDbSnp(snp)
            nummnp += len(subsnps)
            #sys.stderr.write("mnp, %d\n" % len(subsnps))
            for s in subsnps:
                snps.append( s )
        elif snp.isSnp:
            snps.append( snp )
    #sys.stdout.write("TotalSnps\t%d\n" % (len(snps)))
    sys.stderr.write("numMNP: %d\n" %(nummnp))

    return snps

#============ main ==============
dbsnps = readDbSnps(sys.argv[1])
readSnps(sys.argv[2], dbsnps)
