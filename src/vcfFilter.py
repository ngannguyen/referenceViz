#!/usr/bin/env python

#nknguyen soe ucsc edu
#March 29 2012
#
#Filter redundant variants of a sorted Vcf file
#I.e: for all variants with same position (same chrom, same pos), make sure that the reference allele is the same, and collapse all alternative alleles into a non-redundant list
#
import os, sys, re, copy

class Variant():
    def __init__(self, line):
        line = line.strip()
        items = line.split("\t")
        if len(items) < 5:
            raise ValueError("Wrong vcf format. Expected 8 fields, only see %d fields. Line: %s\n" %(len(items), line))
        self.chr = items[0]
        self.pos = int(items[1]) #convert to base 0
        self.id = items[2]
        if self.id == '.':
            self.id = "%s-%d" %(self.chr, self.pos)
        self.ref = items[3]
        sep = ','
        if re.search("/", items[4]):
            sep = "/"
        alts = items[4].rstrip(sep).split(sep)
        #Remove reference allele from observed alleles:
        self.alts = []
        for a in alts:
            if a != self.ref:
                self.alts.append(a)
        
        for i, a in enumerate( self.alts ):
            #if a == "<DEL>" or a == '-':
            if a == '-':
                self.alts[i] = '<DEL>'

    def setId(self, id):
        self.id = id

    def __cmp__(self, other):
        if self.chr.lstrip('chr') != other.chr.lstrip('chr'):
            return cmp(self.chr.lstrip('chr'), other.chr.lstrip('chr'))
        elif self.pos != other.pos:
            return cmp(self.pos, other.pos)
        else:
            return cmp(self.ref, other.ref)

    def getStr(self):
        return "%s\t%d\t%s\t%s\t%s" %(self.chr, self.pos, self.id, self.ref, ','.join(self.alts))

def collapseVariants(variants):
    if len(variants) == 1:
        return variants[0]
    elif len(variants) == 0:
        return None
    
    variant = copy.copy( variants[0] )
    #Make sure that all the variants are 'equal':
    for i in xrange( 1, len(variants) ):
        if variants[i] != variants[i - 1] or variants[i].ref != variants[i-1].ref:
            raise ValueError("Trying to collapse non-redundant variants: %s and %s\n" %( variants[i].getStr(), variants[i-1].getStr() ) )
    
    #Now collapse:
    #id = "%s-%d" %(variant.chr, variant.pos)
    for i in xrange( 1, len(variants) ):
        for a in variants[i].alts:
            if a not in variant.alts:
                variant.alts.append(a)
        #if variants[i].id != id:
        #    variant.id += ",%s" %variants[i].id
    return variant

def getRSlist(id):
    snps = []
    if re.search("dbsnp", id):
        items = id.split('-')
        for item in items:
            if re.search("dbsnp", item):
                items2 = item.split('_')
                if len(items2) >= 2 and re.search("rs", items2[1]):
                    snps.append(items2[1])
    else:
        snps = [id]
    return snps

def readVcf(f):
    visitedVariants = {}
    variants = []
    prev = None
    currVars = [] #list of current redundant variants
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        v = Variant(line)

        #Ignore snps rs***** that have already beed added to the variants list
        if re.search("rs", v.id):
            ids = getRSlist(v.id)
            visited = False
            for id in ids:
                if id in visitedVariants: 
                    visitedVariants[id] += 1
                    visited = True
                else:
                    visitedVariants[id] = 1
            if visited:
                continue

        if not prev or prev == v:
            currVars.append( v )
        else:
            #Collapse currVars:
            collapsedVariant =  collapseVariants( currVars )
            #collapsedVariant.setId( str(len(variants)) )
            variants.append(collapsedVariant)
            #Reset currVars
            currVars = [v]
        prev = v
    if len(currVars) > 0:
        collapsedVariant =  collapseVariants( currVars )
        collapsedVariant.setId( str(len(variants)) )
        variants.append(collapsedVariant)
        #variants.append( collapseVariants(currVars) )
    sys.stderr.write( "Total visited variants: %d; redundant: %d\n" %(len(visitedVariants), sum(visitedVariants.values())) )
    return variants

def writeVcf(f, variants):
    for v in variants:
        f.write("%s\n" %v.getStr())

#====== main =====
variants = readVcf(sys.stdin)
writeVcf(sys.stdout, variants)





