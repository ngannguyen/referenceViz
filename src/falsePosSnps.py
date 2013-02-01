#!/usr/bin/env python

"""
nknguyen soe ucsc edu
May 22 2012
This script examines potential false positive SNPs in the n->1 problem, where duplications exist in the population but not the reference, and reads of the absent copy falsely mapped to the copy the reference has.

Particularly, we will use trio data and linkage to detect if a SNP is more likely to be false or not
Input: a/ file that contains regions of interest
       b/ phased trio SNPs calls (from 1k GP), which should include SNPs calls for the father, mother and child
Output: 
"""

import os, re, sys


class TrioSnp():
    def __init__(self, line, ref, alts):
        alleleStr = line.split(':')[0]
        
        self.phased = False
        items = alleleStr.split('/')
        if re.search("\|", alleleStr):
            self.phased = True
            items = alleleStr.split('|')
        if len(items) != 2:
            raise ValueError("SNP does not have two alleles: %s" % line)
        
        allAlleles = [ref]
        allAlleles.extend(alts)
        self.alleles = [ allAlleles[int(i)] for i in items ]
        self.het = True
        if self.alleles[0] == self.alleles[1]:
            self.het = False

class Snp():
    def __init__(self, line):
        line = line.strip()
        items = line.split('\t')
        self.chr = items[0].lstrip('chr')
        self.pos = int(items[1])
        self.id = items[2]
        self.ref = items[3]
        self.alts = items[4].split(',')
        self.qual = items[5]
        self.filter = items[6]
        self.father = TrioSnp(items[9], self.ref, self.alts)
        self.mother = TrioSnp(items[10], self.ref, self.alts)
        self.child = TrioSnp(items[11], self.ref, self.alts)

def getPhasedSnps(pos2snp):
    pos2phasedSnp = {}
    for pos, snp in pos2snp.iteritems():
        if snp.father.phased and snp.mother.phased and snp.child.phased:
            pos2phasedSnp[pos] = snp
    return pos2phasedSnp

def checkUnphased(snp, k, total, totalsnps, selected, selectedsnps, chr2pos):
    if k == 'f':
        trioSnp = snp.father
    elif k == 'm':
        trioSnp = snp.mother
    elif k == 'c':
        trioSnp = snp.child
    else:
        raise KeyError('Key %s is unknown' %k)
    if trioSnp.alleles[0] == snp.ref and trioSnp.alleles[1] == snp.ref:
        return
    totalsnps[k] += 1
    if snp.chr in chr2pos and snp.pos in chr2pos[snp.chr]:
        selectedsnps[k] += 1
    if not trioSnp.phased:
        total[k] += 1
        if snp.chr in chr2pos and snp.pos in chr2pos[snp.chr]:
            selected[k] += 1

def getRatio(numer, denom):
    ratio = {}
    for k, d in denom.iteritems():
        if k not in numer:
            continue
        if d == 0:
            ratio[k] = 0.0
        else:
            ratio[k] = numer[k]*100.0/d
    return ratio

def countHaps(chr2pos, chr2pos2snp):
    numhap2count = {}
    selectedNumhap2count = {}
    
    for chr in chr2pos2snp:
        pos2snp = chr2pos2snp[chr]
        pos2phasedSnps = getPhasedSnps(pos2snp) #phased SNPs = SNPs where both parents & child are phased
        sortedPos = sorted(pos2phasedSnps.keys())
        pos2index = {}
        for i, p in enumerate(sortedPos):
            pos2index[p] = i

        for pos, snp in pos2snp.iteritems():
            if pos not in pos2phasedSnps:
                continue
            i = pos2index[pos]
            if i == 0 or i == len(sortedPos) - 1:
                continue

            lefti = i - 1
            leftpos = sortedPos[lefti]
            #while leftpos in chr2pos[chr].keys() and pos - leftpos <= 5000 and lefti > 0:
            #    lefti -= 1
            #    leftpos = sortedPos[lefti]
            #if leftpos in chr2pos[chr].keys() or pos - leftpos > 5000:
            #    continue
            
            righti = i + 1
            rightpos = sortedPos[righti]
            #while rightpos in chr2pos[chr].keys() and rightpos - pos <= 5000 and righti < len(sortedPos) - 1:
            #    righti += 1
            #    rightpos = sortedPos[righti]
            #if rightpos in chr2pos[chr].keys() or rightpos - pos > 5000:
            #    continue

            left = pos2phasedSnps[ sortedPos[lefti] ]
            right = pos2phasedSnps[ sortedPos[righti] ]
           
            #HACK
            if not left.father.het or not left.mother.het or not left.child.het or \
               not snp.father.het or not snp.mother.het or not snp.child.het or \
               not right.father.het or not right.mother.het or not right.child.het:
               continue
            #END HACK

            haps = []
            for i in [0, 1]:
                hapi = [left.mother.alleles[i], snp.mother.alleles[i], right.mother.alleles[i]]
                if hapi not in haps:
                    haps.append(hapi)
            
            for i in [0, 1]:
                hapi = [left.father.alleles[i], snp.father.alleles[i], right.father.alleles[i]]
                if hapi not in haps:
                    haps.append(hapi)

            numhap = len(haps)
            if numhap not in numhap2count:
                numhap2count[ numhap ] = 1
            else:
                numhap2count[ numhap ] += 1

            if chr in chr2pos and pos in chr2pos[chr]:
                if numhap not in selectedNumhap2count:
                    selectedNumhap2count[ numhap ] = 1
                else:
                    selectedNumhap2count[ numhap ] += 1

    #Print summary:
    total = sum([ count for count in numhap2count.values() ])
    selected = sum([ count for count in selectedNumhap2count.values() ])
    for num, count in numhap2count.iteritems():
        numhap2count[num] = numhap2count[num]*100.0/total
    for num, count in selectedNumhap2count.iteritems():
        selectedNumhap2count[num] = selectedNumhap2count[num]*100.0/selected
    cats = sorted(numhap2count.keys())
    print "\t%s" %( '\t'.join([str(c) for c in cats]) )
    sys.stdout.write("Overall")
    for c in cats:
        sys.stdout.write("\t%f" % numhap2count[c])
    sys.stdout.write("\nSelectedRegions")
    for c in cats:
        if c not in selectedNumhap2count:
            sys.stdout.write("\t0.0")
        else:
            sys.stdout.write("\t%f" %selectedNumhap2count[c])

def updatehom2het(snp, k, homAlt, homRef, het, selectedhomAlt, selectedhomRef, selectedhet, chr2pos):
    if k == 'f':
        trioSnp = snp.father
    elif k == 'm':
        trioSnp = snp.mother
    elif k == 'c':
        trioSnp = snp.child
    else:
        raise KeyError('Key %s is unknown' %k)

    #if trioSnp.alleles[0] == snp.ref and trioSnp.alleles[1] == snp.ref:
    #    return
    if trioSnp.het:
        het[k] += 1
        if snp.chr in chr2pos and snp.pos in chr2pos[snp.chr]:
            selectedhet[k] += 1
    elif trioSnp.alleles[0] == snp.ref: #homRef (or AA)
        homRef[k] += 1
        if snp.chr in chr2pos and snp.pos in chr2pos[snp.chr]:
            selectedhomRef[k] += 1
    else:
        homAlt[k] += 1
        if snp.chr in chr2pos and snp.pos in chr2pos[snp.chr]:
            selectedhomAlt[k] += 1

def getfreq(aa, ab, bb):
    total = float(aa + ab + bb)
    return aa/total, ab/total, bb/total

def getallelefreq(aa, ab, bb):
    total = float(aa + ab + bb)
    return (aa + ab/2.0)/total , (bb + ab/2.0)/total

def expectedfreq(a1, b1, a2, b2):
    return a1*a2, a1*b2 + b1*a2, b1*b2

def homo2het(chr2pos, chr2pos2snp):
    homRef = {'f':0, 'm':0, 'c': 0}
    homAlt = {'f':0, 'm':0, 'c': 0}
    het = {'f':0, 'm':0, 'c': 0}
    selectedhomRef = {'f':0, 'm':0, 'c': 0}
    selectedhomAlt = {'f':0, 'm':0, 'c': 0}
    selectedhet = {'f':0, 'm':0, 'c': 0}

    #Total:
    for chr in chr2pos2snp:
        for pos, snp in chr2pos2snp[chr].iteritems():
            updatehom2het(snp, 'f', homAlt, homRef, het, selectedhomAlt, selectedhomRef, selectedhet, chr2pos)
            updatehom2het(snp, 'm', homAlt, homRef, het, selectedhomAlt, selectedhomRef, selectedhet, chr2pos)
            updatehom2het(snp, 'c', homAlt, homRef, het, selectedhomAlt, selectedhomRef, selectedhet, chr2pos)
                            
    #s2t = getRatio(selectedRatio, totalRatio)
    #for k, v in s2t.iteritems():
    #    s2t[k] /= 100.0
    print "Category\tFather\tMother\tChild"
    print "OVERALL"
    print "OverallAB\t%d\t%d\t%d" %(het['f'], het['m'], het['c'])
    print "OverallBB\t%d\t%d\t%d" %(homAlt['f'], homAlt['m'], homAlt['c'])
    print "OverallAA\t%d\t%d\t%d" %(homRef['f'], homRef['m'], homRef['c'])
    
    print "\n\tfreq A\t freq B"
    fA, fB = getallelefreq(homRef['f'], het['f'], homAlt['f'])
    print "Father\t%f\t%f" %(fA, fB)
    mA, mB = getallelefreq(homRef['m'], het['m'], homAlt['m'])
    print "Mother\t%f\t%f" %(mA, mB)
    print "\nChild Freq\tAB\tBB\tAA"
    cAA, cAB, cBB = expectedfreq(fA, fB, mA, mB)
    print "Expected\t%f\t%f\t%f" %(cAB, cBB, cAA)
    ocAA, ocAB, ocBB = getfreq(homRef['c'], het['c'], homAlt['c'])
    print "Observed\t%f\t%f\t%f" %(ocAB, ocBB, ocAA)

    print "\n\nSELECTED"
    print "SelectedAB\t%d\t%d\t%d" %(selectedhet['f'], selectedhet['m'], selectedhet['c'])
    print "SelectedBB\t%d\t%d\t%d" %(selectedhomAlt['f'], selectedhomAlt['m'], selectedhomAlt['c'])
    print "SelectedAA\t%d\t%d\t%d" %(selectedhomRef['f'], selectedhomRef['m'], selectedhomRef['c'])
    
    print "\n\tfreq A\t freq B"
    fA, fB = getallelefreq(selectedhomRef['f'], selectedhet['f'], selectedhomAlt['f'])
    print "Father\t%f\t%f" %(fA, fB)
    mA, mB = getallelefreq(selectedhomRef['m'], selectedhet['m'], selectedhomAlt['m'])
    print "Mother\t%f\t%f" %(mA, mB)
    print "\nChild Freq\tAB\tBB\tAA"
    cAA, cAB, cBB = expectedfreq(fA, fB, mA, mB)
    print "Expected\t%f\t%f\t%f" %(cAB, cBB, cAA)
    ocAA, ocAB, ocBB = getfreq(selectedhomRef['c'], selectedhet['c'], selectedhomAlt['c'])
    print "Observed\t%f\t%f\t%f" %(ocAB, ocBB, ocAA)

    
    #totalRatio = getRatio(homo, het)
    #selectedRatio = getRatio(selectedhomo, selectedhet)
    #print "Overall\t%f\t%f\t%f"  %(totalRatio['f']/100.0, totalRatio['m']/100.0, totalRatio['c']/100.0)
    #print "Selected\t%f\t%f\t%f"  %(selectedRatio['f']/100.0, selectedRatio['m']/100.0, selectedRatio['c']/100.0)
    #print "S/T\t%f\t%f\t%f"  %(s2t['f'], s2t['m'], s2t['c'])


def checkUnphasedSnps(chr2pos, chr2pos2snp):
    total = {'f':0, 'm':0, 'c': 0}
    totalsnps = {'f':0, 'm':0, 'c': 0}
    selected = {'f':0, 'm':0, 'c': 0}
    selectedsnps = {'f':0, 'm':0, 'c': 0}

    #Total:
    for chr in chr2pos2snp:
        for pos, snp in chr2pos2snp[chr].iteritems():
            checkUnphased(snp, 'f', total, totalsnps, selected, selectedsnps, chr2pos)
            checkUnphased(snp, 'm', total, totalsnps, selected, selectedsnps, chr2pos)
            checkUnphased(snp, 'c', total, totalsnps, selected, selectedsnps, chr2pos)
                            
    print "Unphased SNPs: Father\tMother\tChild"
    totalRatio = getRatio(total, totalsnps)
    selectedRatio = getRatio(selected, selectedsnps)
    s2t = getRatio(selectedRatio, totalRatio)
    for k, v in s2t.iteritems():
        s2t[k] /= 100.0
    print "Total\t%f\t%f\t%f"  %(totalRatio['f'], totalRatio['m'], totalRatio['c'])
    print "Selected\t%f\t%f\t%f"  %(selectedRatio['f'], selectedRatio['m'], selectedRatio['c'])
    print "S/T\t%f\t%f\t%f"  %(s2t['f'], s2t['m'], s2t['c'])
    

def checkTrioSnps(chr2pos, chr2pos2snp):
    #chr2pos is a list of suspicious positions
    #chr2pos2snp is the list of 1k genome snps of the trios (father, mother and child)
    #for each suspicious position, check to see if there is a phased snp called by 1kGP
    #If yes, look for the closest left phased snp, and closest right phased snp
    #and check to see if the child haplotype requires any recombination, if yes and the 
    #distance between left, current, and right SNPs are very close, it is probably a spurious SNP
    
    numRecoms = 0 #number of cases where recombination is needed to explain 
    numExams = 0 #number of snp tested
    range = 1000
    for chr in chr2pos:
        if chr not in chr2pos2snp:
            continue
        pos2snp = chr2pos2snp[chr]
        pos2phasedSnps = getPhasedSnps(pos2snp) #phased SNPs = SNPs where both parents & child are phased
        sortedPos = sorted(pos2phasedSnps.keys())
        pos2index = {}
        for i, p in enumerate(sortedPos):
            pos2index[p] = i

        for pos in chr2pos[chr].keys():
            if pos not in pos2phasedSnps:
                continue
            snp = pos2snp[pos]
            if snp.child.alleles[0] == snp.ref and snp.child.alleles[1] == snp.ref: #no SNP in the child
                continue

            i = pos2index[pos]
            if i < 2 or i > len(pos2index) - 2:
                continue

            lefti = i - 1
            leftpos = sortedPos[lefti]
            llpos = sortedPos[lefti - 1]
            #while leftpos in chr2pos[chr].keys() and pos - leftpos <= 1000 and lefti > 0:
            #    lefti -= 1
            #    leftpos = sortedPos[lefti]
            #if leftpos in chr2pos[chr].keys() or pos - leftpos > 1000:
            #    continue

            righti = i + 1
            rightpos = sortedPos[righti]
            rrpos = sortedPos[righti + 1]
            #while rightpos in chr2pos[chr].keys() and rightpos - pos <= 1000 and righti < len(sortedPos) - 1:
            #    righti += 1
            #    rightpos = sortedPos[righti]
            #if rightpos in chr2pos[chr].keys() or rightpos - pos > 1000:
            #    continue
            
            left = pos2phasedSnps[ leftpos ]
            right = pos2phasedSnps[ rightpos ]
            ll = pos2phasedSnps[ llpos ]
            rr = pos2phasedSnps[ rrpos ]
            numExams += 1
            
            childhaps = [] 
            for i in [0, 1]:
                #hapi = [left.child.alleles[i], snp.child.alleles[i], right.child.alleles[i]]
                hapi = [ll.child.alleles[i], left.child.alleles[i], snp.child.alleles[i], right.child.alleles[i], rr.child.alleles[i]]
                childhaps.append(hapi)
            
            motherhaps = [] 
            for i in [0, 1]:
                #hapi = [left.mother.alleles[i], snp.mother.alleles[i], right.mother.alleles[i]]
                hapi = [ll.mother.alleles[i], left.mother.alleles[i], snp.mother.alleles[i], right.mother.alleles[i], rr.mother.alleles[i]]
                motherhaps.append(hapi)
            
            fatherhaps = [] 
            for i in [0, 1]:
                #hapi = [left.father.alleles[i], snp.father.alleles[i], right.father.alleles[i]]
                hapi = [ll.father.alleles[i], left.father.alleles[i], snp.father.alleles[i], right.father.alleles[i], rr.father.alleles[i]]
                fatherhaps.append(hapi)
            
            print ll.pos, left.pos, snp.pos, right.pos, rr.pos
            print motherhaps
            print fatherhaps
            print childhaps

            ok = False
            if (childhaps[0] in motherhaps and childhaps[1] in fatherhaps) or (childhaps[0] in fatherhaps and childhaps[1] in motherhaps):
                ok = True
            if not ok:
                numRecoms += 1
    sys.stdout.write("Number of SNPs examed: %d\n" %numExams)
    sys.stdout.write("Potential spurious SNPs: %d\t%f%%\n" %(numRecoms, numRecoms*100.0/numExams))

def readPosFile(file):
    f = open(file, 'r')
    chr2pos = {}
    for line in f:
        items = line.split('\t')
        chritems = items[0].split('.')
        #hg18.chr6.170899992.28585733.4970599.1
        chr = chritems[1].lstrip('chr')
        offset = int(chritems[3])
        pos = int(items[1]) + offset
        if chr not in chr2pos:
            chr2pos[chr] = {pos: 1}
        else:
            chr2pos[chr][pos] = 1
    f.close()
    return chr2pos

def readVcfFile(file):
    chr2pos2snp = {}
    f = open(file, 'r')
    for line in f:
        if line[0] == '#':
            continue
        snp = Snp(line)
        if snp.chr not in chr2pos2snp:
            chr2pos2snp[snp.chr] = { snp.pos: snp }
        else:
            chr2pos2snp[snp.chr][snp.pos] = snp
    f.close()
    return chr2pos2snp 

def main():
    vcffile = sys.argv[1]
    posfile = sys.argv[2]
    chr2pos = readPosFile(posfile)
    chr2pos2snp = readVcfFile(vcffile)
    #checkTrioSnps(chr2pos, chr2pos2snp)
    #checkUnphasedSnps(chr2pos, chr2pos2snp)
    #countHaps(chr2pos, chr2pos2snp)
    homo2het(chr2pos, chr2pos2snp)


if __name__ == '__main__':
    main()


