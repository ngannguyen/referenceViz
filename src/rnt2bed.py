#!/usr/bin/env python

'''
Thu Jul 18 10:06:47 PDT 2013
Convert files in rnt to bed format. Also change the names to alpha-numeric only
Indir/
    Sample_1/
        chr_1.rnt
        ...
    Sample_2/
        ...
    ...
Outdir/
    Sample1/
        chr1.bed
'''

import os, sys, re, time
from sonLib.bioio import system
from optparse import OptionParser

def convertToAlnum(inputStr):
    items = re.split("[^a-zA-Z0-9]", inputStr)
    items = [item.capitalize() for item in items] 
    newStr = "".join(items)
    if not newStr.isalnum():
        print newStr
    return newStr

def convert2bedfile(infile, outfile, chr, color):
    #Escherichia coli 042, complete genome - 1..5241977
    #114 RNAs
    #Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
    #228620..230161  +       1542    387605479       -       -       -       -       16S ribosomal RNA
    #230230..230303  +       74      387605479       -       -       -       -       Anticodon: GAT
    ifh = open(infile, 'r')
    ofh = open(outfile, 'w')
    ifh.readline()
    ifh.readline()
    ifh.readline()
    for line in ifh:
        items = line.strip().split('\t')
        assert(len(items)) == 9
        locations = items[0].split('..')
        start = int(locations[0]) -1 #convert to base 0
        end = int(locations[1]) #end is exclusive
        strand = items[1]
        #l = int(items[2])
        pid = items[3]
        gene = items[4]
        synonym = items[5]
        product = items[8]
        if gene != '-':
            name = gene
        elif product != '-':
            productItems = re.split("[^a-zA-Z0-9]+", product)
            name = '-'.join(productItems)
        else:
            name = synonym
        #bed: chrom chromstart chromend name score strand thickstart thickend itemRgb
        ofh.write("%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t%s\n" %(chr, start, end, name, strand, start, end, color))

    ifh.close()
    ofh.close()

def convert2beds(indir, outdir, color, genome2seqs, names):
    old2new = {}
    samples = os.listdir(indir)
    if not os.path.exists(outdir):
        system("mkdir -p %s" %outdir)

    for sample in samples:
        sampledir = os.path.join(indir, sample)
        #convert sample to alnum only:
        outsample = convertToAlnum(sample)
        seqs = []
        if outsample.lower() in genome2seqs:
            seqs = genome2seqs[outsample.lower()]
            outsample = names[outsample.lower()]
        outsampledir = os.path.join(outdir, outsample)
        system("mkdir -p %s" %outsampledir)
        if outsample != sample:
            old2new[sample] = outsample

        chrfiles = os.listdir(sampledir)
        for chrfile in chrfiles:
            chritems = chrfile.split('.')
            if len(chritems) == 1:
                chr = chrfile
            else:
                chr = '.'.join(chritems[:-1])
            
            outchr = convertToAlnum(chr)
            for s in seqs:
                if re.search(outchr, s):
                    outchr = s

            if outchr != chr:
                old2new[chr] = outchr

            inchrfile = os.path.join(sampledir, chrfile)
            outchrfile = os.path.join(outsampledir, "%s.bed" %outchr)
            convert2bedfile(inchrfile, outchrfile, outchr, color)
    
    for o, n in old2new.iteritems():
        sys.stdout.write("%s\t%s\n" %(o, n))

def readGenomeSeqs(file):
    genome2seqs = {}
    lcname2name = {}

    f = open(file, 'r')
    genome = ''
    seqs = []
    for l in f:
        l = l.strip()
        if len(l) == 0:
            continue
        if l[0] == '#':
            if genome != '' and len(seqs) > 0:
                genome2seqs[genome.lower()] = seqs
                lcname2name[genome.lower()] = genome
            genome = l.split()[-1]
            seqs = []
        else:
            seqs.extend(l.split(','))
            
    if genome != '' and len(seqs) > 0:
        genome2seqs[genome.lower()] = seqs
        lcname2name[genome.lower()] = genome
    f.close()
    return genome2seqs, lcname2name

def main():
    usage = '%prog [options] indir outdir'
    parser = OptionParser(usage=usage)
    parser.add_option('--color', dest='color', default='0,150,0', help='itemRGB to color the bed tracks. Default=%default')
    parser.add_option('--genomeSequences', dest='genomeSeqs', help='If specified, try to convert the sequence names to match with names in this file')
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Require two input arguments <indir> and <outdir>. Only see %d\n" %len(args))
    indir = args[0]
    outdir = args[1]
    genome2seqs = {}
    names = {}
    if options.genomeSeqs:
        genome2seqs, names = readGenomeSeqs(options.genomeSeqs)
    
    convert2beds(indir, outdir, options.color, genome2seqs, names)    

if __name__ == '__main__':
    main()

