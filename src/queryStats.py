#!/usr/bin/env python

"""
"""
import os, sys, re

#Read transcriptIDtoGeneName
def mapTx2gene(file):
    gene2tx = {}
    tx2gene = {}
    f = open(file, 'r')
    for line in f:
        if line[0] == '#':
            continue
        items = line.strip().split('\t')
        if len(items) != 2:
            continue
        #Update gene2tx:
        if items[1] not in gene2tx:
            gene2tx[items[1]] = [items[0]]
        elif items[0] not in gene2tx[items[1]]:
            gene2tx[items[1]].append(items[0])
        
        #Update tx2gene:
        if items[0] not in tx2gene:
            tx2gene[items[0]] = items[1]
    f.close()
    return gene2tx, tx2gene

#Read queryStatsFile:
def readQueryStats(file, tx2gene):
    gene2tx2cov = {} #key = gene, val = {transcriptID: maxQcover}
    f = open(file, 'r')
    for line in f:
        if line[0] == '#':
            continue
        items = line.strip().split('\t')
        if len(items) <= 10:
            continue
        txid = items[0].split('.')[0]
        maxQCover = float(items[7])
        if txid not in tx2gene:
            sys.stderr.write("Transcript without any gene matched: %s\n" %txid)
            continue
        gene = tx2gene[txid]
        if gene not in gene2tx2cov:
            gene2tx2cov[gene] = {txid: maxQCover}
        elif txid not in gene2tx2cov[gene] or maxQCover > gene2tx2cov[gene][txid]:
            gene2tx2cov[gene][txid] = maxQCover

    f.close()
    return gene2tx2cov

#Counting:
def getCounts(gene2tx2cov, gene2tx):
    totalGenes = len(gene2tx2cov.keys())
    stats = {'noBrokenTx':0, 'atLeastOneIntactTx':0}
    genes = {'noBrokenTx':[], 'atLeastOneIntactTx':[]}

    for g in gene2tx2cov:
        totaltx = len( gene2tx[g] )
        maptx = len( gene2tx2cov[g] )
        if totaltx != maptx:
            sys.stderr.write("Gene with inconsistent number of transcripts: %s, totalFromTx2gene: %d, totalfromPsl: %d\n" %(g, totaltx, maptx))
        numfulltx = 0
        for tx in gene2tx2cov[g]:
            #if gene2tx2cov[g][tx] >= 1.0:
            if gene2tx2cov[g][tx] >= 0.95:
                numfulltx += 1
        if numfulltx == maptx:
            stats['noBrokenTx'] += 1
            genes['noBrokenTx'].append(g)
        elif numfulltx >= 1:
            stats['atLeastOneIntactTx'] += 1
            genes['atLeastOneIntactTx'].append(g)

    sys.stdout.write("Total Genes\t%d\n" %(totalGenes))
    sys.stdout.write("No Broken tx\t%d\n" %(stats['noBrokenTx']))
    sys.stdout.write("At least one unbroken tx\t%d\n" %(stats['atLeastOneIntactTx']))
    sys.stdout.write("Others\t%d\n" %( totalGenes - stats['atLeastOneIntactTx'] - stats['noBrokenTx']))

    #
    for g in gene2tx2cov:
        if g not in genes['noBrokenTx'] and g not in genes['atLeastOneIntactTx']:
            sys.stderr.write("%s\n"%g)
    return genes

#main:
tx2genefile = sys.argv[1]
crefstats = sys.argv[2]
hg19stats = sys.argv[3]

gene2tx, tx2gene = mapTx2gene(tx2genefile)
crefGene2tx2cov = readQueryStats(crefstats, tx2gene)
sys.stdout.write("\nCref\n")
crefGenes = getCounts(crefGene2tx2cov, gene2tx)

hg19Gene2tx2cov = readQueryStats(hg19stats, tx2gene)
sys.stdout.write("\nhg19\n")
hg19Genes = getCounts(hg19Gene2tx2cov, gene2tx)

sys.stderr.write("\nAll transcripts fully mapped to CREF, None to hg19:\n")
for g in crefGenes['noBrokenTx']:
    if g not in hg19Genes['noBrokenTx'] and g not in hg19Genes['atLeastOneIntactTx']:
        sys.stderr.write("%s\n" %g)
    
sys.stderr.write("\nAt least one tx fully mapped to CREF, None to hg19:\n")
for g in crefGenes['atLeastOneIntactTx']:
    if g not in hg19Genes['noBrokenTx'] and g not in hg19Genes['atLeastOneIntactTx']:
        sys.stderr.write("%s\n" %g)

sys.stderr.write("\n\nAll transcripts fully mapped to hg19, None to CREF:\n")
for g in hg19Genes['noBrokenTx']:
    if g not in crefGenes['noBrokenTx'] and g not in crefGenes['atLeastOneIntactTx']:
        sys.stderr.write("%s\n" %g)

sys.stderr.write("\nAt least one tx fully mapped to hg19, None to CREF:\n")
for g in hg19Genes['atLeastOneIntactTx']:
    if g not in crefGenes['noBrokenTx'] and g not in crefGenes['atLeastOneIntactTx']:
        sys.stderr.write("%s\n" %g)








