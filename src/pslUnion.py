#!/usr/bin/env python

'''
Fri Feb 22 13:05:44 PST 2013
nknguyen soe ucsc edu
Intersect multiple lists of genes (in psl formats)

Input:
   Directory contains input fasta files, one for each sample
   output directory
   sample.lst  : list of samples(/species/individuals/...) involved

Output:
   File containing the union genelist in the format:
<geneID>   <comma,sep,list,of,homologous,genes>
where geneID is assigned by the program uniquely for each gene family.
'''

import copy, os, sys, re, time, random, gzip
from optparse import OptionParser
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import system

######################### Obj classes #####################
class Psl():
    '''Psl record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) != 21: 
            raise PslFormatError("Psl format requires 21 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        self.desc = line
        self.matches = int(items[0])
        self.misMatches = int(items[1])
        self.repMatches = int(items[2])
        self.nCount = int(items[3])
        self.qNumInsert = int(items[4])
        self.qBaseInsert = int(items[5]) #number of bases inserted in query
        self.tNumInsert = int(items[6]) #number of inserts in target
        self.tBaseInsert = int(items[7]) #number of bases inserted in target
        self.strand = items[8] #query strand
        self.qName = items[9] 
        self.qSize = int(items[10])
        self.qStart = int(items[11]) #base 0
        self.qEnd = int(items[12])
        self.tName = items[13]
        self.tSize = int(items[14])
        self.tStart = int(items[15])
        self.tEnd = int(items[16])
        self.blockCount = int(items[17])
        self.blockSizes = [int(s) for s in items[18].rstrip(',').split(',')]
        self.qStarts = [int(s) for s in items[19].rstrip(',').split(',')]
        self.tStarts = [int(s) for s in items[20].rstrip(',').split(',')]
        if len(self.blockSizes) != self.blockCount or len(self.qStarts) != self.blockCount or len(self.tStarts) != self.blockCount:
            raise PslFormatError("Psl format requires that the number of items in blockSizes, qStarts, tStarts is equal to blockCount. Line: %s\n" %line)

class GeneFamily():
    '''Representing a gene family
    '''
    def __init__(self, id):
        self.id = id
        self.genes = []
        self.sample2genes = {}
        self.length = 0
        self.longestGene = None

    def addGene(self, gene, length):
        if gene not in self.genes:
            self.genes.append(gene)
            if length > self.length:
                self.length = length
                self.longestGene = gene

    def addGenes(self, genes, lengths):
        for i, g in enumerate(genes):
            self.addGene(g, lengths[i])

    def addGene2(self, gene, length, sample):
        self.addGene(gene, length)
        if sample not in self.sample2genes:
            self.sample2genes[sample] = [gene]
        elif gene not in self.sample2genes[sample]:
            self.sample2genes[sample].append(gene)

    def addGenes2(self, genes, lengths, sample):
        for i, g in enumerate(genes):
            self.addGene2(g, lengths[i], sample)

    def addGene3(self, gene, sample):
        if gene not in self.genes:
            self.genes.append(gene)
        if sample not in self.sample2genes:
            self.sample2genes[sample] = [gene]
        elif gene not in self.sample2genes[sample]:
            self.sample2genes[sample].append(gene)

    def addGenes3(self, genes, sample):
        for i, g in enumerate(genes):
            self.addGene3(g, sample)
        
    def setId(self, id):
        self.id = id

############################ PIPELINE ############################
class Setup(Target):
    def __init__(self, fadir, samples, outfile, options):
        Target.__init__(self)
        self.fadir = fadir
        self.samples = samples
        self.outfile = outfile
        self.options = options

    def run(self):
        globalTempDir = self.getGlobalTempDir()
        unionDir = os.path.join(globalTempDir, "selfAlign")
        system("mkdir -p %s" %unionDir)
        
        unionFaDir = os.path.join(unionDir, "fa")
        system("mkdir -p %s" %unionFaDir)
        
        unionPickleDir = os.path.join(unionDir, "pickle")
        system("mkdir -p %s" %unionPickleDir)


        #Self-align and Collapse genes into gene families for each sample:
        for sample in self.samples:
            sampleFadir = os.path.join(self.fadir, sample)
            self.addChildTarget( SelfAlign(sampleFadir, unionDir, self.options.prot, self.options.identity, self.options.coverage) )
        
        #Got: globalTempDir/selfAlign/fa containing the non-redundant gene family sequences
        # and globalTempDir/selfAlign/pickle containing the gene family list
        #Now merge samples and get the union of the gene family lists
        self.setFollowOnTarget( UnionGeneLists(None, unionPickleDir, unionFaDir, self.outfile, self.options) )
        
class SelfAlign(Target):
    '''Input: fasta directory containning fasta files of a sample. Self align this sample with itself.
       Collapse all "homologous" sequences (genes) into families.
       Return fasta file of the non-redundant sequences
    '''
    def __init__(self, fadir, outdir, prot, identity, coverage):
        Target.__init__(self)
        self.fadir = fadir
        self.outdir = outdir
        self.prot = prot
        self.identity = identity
        self.coverage = coverage

    def run(self):
        localDir = self.getLocalTempDir()

        #Blat the sequences and return output to localdir/psls
        pslDir = os.path.join(localDir, "psls")
        system("mkdir %s" %pslDir)
        seqfiles = [os.path.join(self.fadir, file)  for file in os.listdir(self.fadir) ]
        for file1 in seqfiles:
            for file2 in seqfiles:
                outfile = os.path.join(pslDir, "%s-%s.psl" %(os.path.basename(file1), os.path.basename(file2)))
                #cmd = "blat %s %s %s -minIdentity=90 -noHead " %(file1, file2, outfile)
                cmd = "blat %s %s %s -minIdentity=%f -noHead " %(file1, file2, outfile, self.identity)
                if self.prot:
                    cmd += " -prot -minScore=10 "
                system(cmd)
        
        #Aggregate & filter all the alignments:
        sample = os.path.basename(self.fadir)
        pslfile = os.path.join(self.getGlobalTempDir(), "%s.psl" % sample)
        #cmd = "cat %s/*psl | awk '$1/$11 >= 0.90 {print $0}' > %s" %(pslDir, pslfile)
        cmd = "cat %s/*psl | awk '$1/$11*100.0 >= %f && $1/$15*100.0 >= %f {print $0}' > %s" %(pslDir, self.coverage, self.coverage, pslfile)
        system(cmd)

        #CollapseDupGenes:
        self.addChildTarget( CollapseDupGenes(sample, pslfile, self.outdir, seqfiles) )

class CollapseDupGenes(Target):
    '''Read in self-alignment psl file of a sample, cluster homologous genes into gene family
       Return fasta file of the non-redundant sequences and pickle file of the union fams
    '''
    def __init__(self, sample, infile, outdir, fafiles):
        Target.__init__(self)
        self.sample = sample
        self.infile = infile
        self.outdir = outdir
        self.fafiles = fafiles

    def run(self):
        name2psls = readPslFile(self.infile, False)
        fams = collapseDups(name2psls, self.sample)
       
        #Extract the non-redundant (/union) fasta
        fadir = os.path.join(self.outdir, "fa")
        unionFaFile = os.path.join(fadir, "%s.fa" %self.sample)
        header2seq = readFastas(self.fafiles)
        writeFasta2(header2seq, fams, unionFaFile) 
        #system("cp %s %s" %(unionFaFile, os.getcwd()))

        #Pickle the fams
        unionPickleDir = os.path.join(self.outdir, "pickle")
        pickleFile = os.path.join(unionPickleDir, "%s.pickle" %self.sample)
        pickle.dump( fams, gzip.open(pickleFile, "wb") )

class UnionGeneLists(Target):
    '''Divide and conquer, set up jobs to merge 2 gene lists at a time
    '''
    def __init__(self, prevPickledir, pickledir, fadir, outfile, options):
        Target.__init__(self)
        self.prevPickledir = prevPickledir
        self.pickledir = pickledir #directory containning pickle files of current gene family lists
        self.fadir = fadir
        self.outfile = outfile
        self.options = options

    def run(self):
        #if self.prevPickledir:
        #    system("rm -R %s" %self.prevPickledir)
        
        inPickles = os.listdir(self.pickledir)
        numSam = len(inPickles)

        if numSam == 1: #done merging, print to output file
            fams = pickle.load(gzip.open(os.path.join(self.pickledir, inPickles[0]), "rb"))
            printUnionList(fams, self.outfile)
            panNcoreGenes(fams, "%s-panNcore" %self.outfile, self.options.samples)
        else: 
            globalTempDir = self.getGlobalTempDir()
            pickleOutdir = os.path.join(globalTempDir, "%d" % numSam, "pickle")
            system("mkdir -p %s" % pickleOutdir)
            
            faOutdir = os.path.join(globalTempDir, "%d" % numSam, "fa")
            system("mkdir -p %s" % faOutdir)
            
            if numSam % 2 == 1: #odd number of files
                oddSample = inPickles.pop().rstrip("pickle").rstrip(".")
                system("mv %s %s" %(os.path.join(self.pickledir, "%s.pickle" % oddSample), pickleOutdir))
                system("mv %s %s" %(os.path.join(self.fadir, "%s.fa" % oddSample), faOutdir))
            assert len(inPickles) %2 == 0
            
            for i in xrange(0, len(inPickles), 2):
                sample1 = inPickles[i].rstrip("pickle").rstrip(".")
                sample2 = inPickles[i+1].rstrip("pickle").rstrip(".")
                pickle1 = os.path.join(self.pickledir, inPickles[i])
                pickle2 = os.path.join(self.pickledir, inPickles[i+1])
                
                fa1 = os.path.join(self.fadir, "%s.fa" % sample1 ) 
                fa2 = os.path.join(self.fadir, "%s.fa" % sample2)
                assert os.path.exists(fa1)
                assert os.path.exists(fa2)
                
                #outPickleFile = os.path.join( pickleOutdir, "%d-%d-%s-%s.pickle" %(numSam, i, sample1, sample2) )
                #outFaFile = os.path.join( faOutdir, "%d-%d-%s-%s.fa" %(numSam, i, sample1, sample2) )
                outPickleFile = os.path.join( pickleOutdir, "%d-%d.pickle" %(numSam, i) )
                outFaFile = os.path.join( faOutdir, "%d-%d.fa" %(numSam, i) )
                self.addChildTarget( UnionPairGeneLists(pickle1, pickle2, fa1, fa2, outPickleFile, outFaFile, self.options) )
                
                #DEBUGGING
                #fams1 = pickle.load(gzip.open(pickle1, "rb"))
                #printUnionList(fams1, "%s_%s" %(self.outfile, os.path.basename(pickle1)))
                #fams2 = pickle.load(gzip.open(pickle2, "rb"))
                #printUnionList(fams2, "%s_%s" %(self.outfile, os.path.basename(pickle2)))
                #END DEBUGGING
            self.setFollowOnTarget( UnionGeneLists(self.pickledir, pickleOutdir, faOutdir, self.outfile, self.options) )

class UnionPairGeneLists(Target):
    '''Find the union of two input gene family lists
    '''
    def __init__(self, pickle1, pickle2, fa1, fa2, pickleOutfile, faOutfile, options):
        Target.__init__(self)
        self.sam1 = os.path.basename(pickle1).split('.')[0]
        self.sam2 = os.path.basename(pickle2).split('.')[0]
        self.name = "%s_%s" %(self.sam1, self.sam2)
        self.fams1 = pickle.load( gzip.open(pickle1, "rb") )
        self.fams2 = pickle.load( gzip.open(pickle2, "rb") )
        self.fa1 = fa1
        self.fa2 = fa2
        self.pickleOutfile = pickleOutfile
        self.faOutfile = faOutfile
        self.options = options

    def run(self):
        #Align fa1 to fa2 and fa2 to fa1
        localDir = self.getLocalTempDir()
        fam2gene2psls2 = pairAlign2(self.fa1, self.fa2, os.path.join(localDir, "%s-%s" %(self.sam1, self.sam2)), self.options.prot, self.options.identity, self.options.coverage)
        fam2gene2psls1 = pairAlign2(self.fa2, self.fa1, os.path.join(localDir, "%s-%s" %(self.sam2, self.sam1)), self.options.prot, self.options.identity, self.options.coverage)
        
        fam2gene2seq1 = readFasta2(self.fa1)
        fam2gene2seq2 = readFasta2(self.fa2)
        unionFam2gene2seq = {}

        #DEBUG
        #ge1 = "gi|209396238|ref|YP_002273063.1|"
        #ge2 = "gi|15804124|ref|NP_290163.1|"

        #Merge the two gene family lists
        unionFams = []
        for i, fam1 in enumerate(self.fams1):
            gene2seq = fam2gene2seq1[fam1.id]
            #logger.info("Fam1 %s; Genes: %s\n" %(fam1.id, ",".join(fam1.genes)))
            if fam1.id not in fam2gene2psls1: #this family does not have a match
                fam = copy.copy(fam1)
                fam.setId( str(len(unionFams)) )
                unionFams.append( fam )
                unionFam2gene2seq[ unionFams[-1].id ] = gene2seq
                #DEBUG
                #if ge1 in gene2seq or ge2 in gene2seq:
                #    logger.info("Fam %s of fams1 does not have any match with any of fams2\n" %(fam1.id))
                #END DEBUG
            else:
                #logger.info( "\tGenesWithHits %s\n" % ",".join(fam2gene2psls1[fam1.id].keys())  )
                merged = False
                for j, fam2 in enumerate(self.fams2):#see if fam1 can merge with any gene family in fams2
                    if fam2.id in fam2gene2psls2 and aligned( fam2gene2psls1[fam1.id], fam2.id, fam2gene2psls2[fam2.id] ):#mergable
                        fam = mergeFamPair(fam1, fam2, str(len(unionFams)))
                        for g2, s2 in fam2gene2seq2[ fam2.id ].iteritems():
                            gene2seq[g2] = s2
                        unionFam2gene2seq[ fam.id ] = gene2seq
                        unionFams.append(fam)
                        self.fams2.remove(fam2)
                        merged = True
                        #DEBUG
                        #if ge1 in gene2seq or ge2 in gene2seq:
                        #    logger.info("\tFam2 %s, Fam1 %s;   Fam2 Genes %s\n" %(fam2.id, fam1.id, ",".join(fam2.genes)))
                        #    logger.info("Fam %s of fams1 (%d) has a reciprocal match with fam %s (%d) of fams2. Merge to Family %s. Len unionFams: %d\n" %(fam1.id, i, fam2.id, j, fam.id, len(unionFams)))
                        #    #logger.info("%d" %fam1.id)
                        #END DEBUG
                        break
                if not merged: #fam1 does not have any reciprocal match with any of fams2, just add it to the union list
                    fam = copy.copy(fam1)
                    fam.setId( str(len(unionFams)) )
                    unionFams.append( fam )
                    unionFam2gene2seq[ unionFams[-1].id ] = gene2seq
                    #DEBUG
                    #if ge1 in gene2seq or ge2 in gene2seq:
                    #    logger.info("Fam %s of fams1 has some match, but NO reciprocal match with any of fams2\n" %(fam1.id))
                    #END DEBUG
                    
        #Whatever left in fams2 do not have reciprocal best match with any in fams1, add it to unionFams:
        for fam2 in self.fams2:
            gene2seq = fam2gene2seq2[fam2.id]
            fam = copy.copy(fam2)
            fam.setId( str(len(unionFams)) )
            unionFams.append( fam )
            unionFam2gene2seq[ unionFams[-1].id ] = gene2seq
            #DEBUG
            #if ge1 in gene2seq or ge2 in gene2seq:
            #    logger.info("Fam %s of fams2 did not have a reciprocal match with fams1\n" %(fam2.id))
            #END DEBUG

        #Pickle unionFams
        pickle.dump(unionFams, gzip.open(self.pickleOutfile, "wb"))
        
        #union fasta file:
        writeFasta3(unionFam2gene2seq, self.faOutfile)

############################ ERROR CLASSES #######################
class PslFormatError(Exception):
    pass

class InputError(Exception):
    pass

class InputFormatError(Exception):
    pass

########################## FUNCTIONS ##############################
def aligned(gene2psls1, fam2id, gene2psls2):
    #ge1 = "gi|209396238|ref|YP_002273063.1|"
    #ge2 = "gi|15804124|ref|NP_290163.1|"
    #if ge1 in gene2psls1 or ge2 in gene2psls1:
    #    logger.info("Aligning Genes %s to Genes %s\n" %( ",".join(gene2psls1.keys()), ",".join(gene2psls2.keys()) ))
    
    #Fam1 can be merged with fam2 if there is at least one gene of fam1 that has match with one gene of fam2
    for gene1, psls1 in gene2psls1.iteritems():#Each gene
        for p1 in psls1:#Each match
            tNameItems = p1.tName.split(';')
            assert len(tNameItems) == 2
            tFam = tNameItems[0]
            tGene = tNameItems[1]
            if tFam == fam2id and tGene in gene2psls2: #this gene of fams2 has some match with fam1
                psls2 = gene2psls2[tGene]
                names1 = [p2.tName for p2 in psls2] #all the matches of gene2 to fams1
                if p1.qName in names1: #reciprocal match
                    return True
    return False

def collapseDups(name2psls, sample):
    fams = []
    gene2fam = {}
    id = 0

    for gene, psls in name2psls.iteritems():
        length = psls[0].qSize
        if gene not in gene2fam:#new family
            id += 1
            fam = GeneFamily( str(id) )
            gene2fam[gene] = fam
            fam.addGene2(gene, length, sample)
            fams.append(fam)
        else:#gene already belongs to a family, now check to see if we can add any of its matches to the family
            fam = gene2fam[gene]

        for p in psls: #for each match
            if p.tName not in gene2fam: #not yet in the family
                tpsls = name2psls[p.tName]
                if gene in [tpsl.tName for tpsl in tpsls]: #is a reciprocal match
                    fam.addGene2(p.tName, p.tSize, sample)
                    gene2fam[p.tName] = fam
    return fams

def pairAlign(fa1, fa2, tempfile, prot, identity, coverage):
    pslfile = "%s-temp" %tempfile
    #cmd = "blat %s %s %s -minIdentity=90 -noHead " %(fa1, fa2, pslfile)
    cmd = "blat %s %s %s -minIdentity=%f -noHead " %(fa1, fa2, pslfile, identity)
    if prot:
        cmd += " -prot -minScore=10 "
    #cmd = "blat %s %s %s -prot -minIdentity=90 -noHead -minScore=10" %(fa1, fa2, pslfile)
    system(cmd)
    
    #Only select matches with >= 90% identity
    #cmd = "awk '$1/$11 >= 0.90 {print $0}' %s > %s" %(pslfile, tempfile)
    cmd = "awk '$1/$11*100.0 >= %f && $1/$15*100.0 >= %f {print $0}' %s > %s" %(coverage, coverage, pslfile, tempfile)
    system(cmd)

    #DEBUG
    #system("cp %s %s" %(tempfile, os.getcwd()))
    #END DEBUG

    name2psls = readPslFile(tempfile, True)
    return name2psls

def pairAlign2(fa1, fa2, tempfile, prot, identity, coverage):
    name2psls = pairAlign(fa1, fa2, tempfile, prot, identity, coverage)
    fam2gene2psls = {}
    for name, psls in name2psls.iteritems():
        items = name.split(';')
        assert len(items) == 2
        fam = items[0]
        gene = items[1]
        if fam not in fam2gene2psls:
            fam2gene2psls[fam] = {gene:psls}
        else:
            fam2gene2psls[fam][gene] = psls
        
        #DEBUG
        #ge1 = "gi|209396238|ref|YP_002273063.1|"
        #ge2 = "gi|15804124|ref|NP_290163.1|"
        #if gene == ge1 or gene == ge2:
        #    logger.info("Gene %s, family %s\n" %(gene, fam))
        #DEBUG

    return fam2gene2psls

def printUnionList(fams, outfile):
    f = open(outfile, 'w')
    f.write("#ID\tGenes\n")
    for fam in fams:
        f.write("%s\t%s\n" %(fam.id, ",".join(fam.genes)))
    f.close()

def mergeFamPair(fam1, fam2, id):
    fam = GeneFamily(id)
    fam.sample2genes = fam1.sample2genes
    fam.genes = fam1.genes

    if fam2.length < fam1.length:
        fam.length = fam1.length
        fam.longestGene = fam1.longestGene
    else:
        fam.length = fam2.length
        fam.longestGene = fam2.longestGene

    for sample, genes in fam2.sample2genes.iteritems():
        fam.addGenes3(genes, sample)
    return fam

def getSample2unionFams(fams):
    sample2fam2genes = {} #key = sample, val = {famid: genes}
    for fam in fams:
        for sample, genes in fam.sample2genes.iteritems():
            if sample not in sample2fam2genes:
                sample2fam2genes[sample] = {fam.id : genes}
            else:
                sample2fam2genes[sample][fam.id] =   genes
    return sample2fam2genes

def panNcoreGenes(fams, outfile, samples):
    '''Add one sample at a time in the order specified by <samples>. 
    Report core gene count and pan gene count each time a new sample added
    '''
    sample2fam2genes = getSample2unionFams(fams)
    
    f = open(outfile, 'w')
    f.write("AddedSample\tSampleCount\tPanGenome\tCoreGenome\tNewGeneFamilies\n")
    panFams = [] #union list of gene families of all samples currently included
    coreFams = [] #core list of gene families shared by all samples currently included

    for i, sample in enumerate(samples):
        newFamCount = 0
        if sample not in sample2fam2genes:
            raise ValueError("Sample %s does not have any gene in the union list!\nSamples in the union list include: %s\n" %(sample, ",".join(sample2fam2genes.keys()) ))
        fam2genes = sample2fam2genes[sample]
        if i == 0:
            coreFams = fam2genes.keys()
            panFams = fam2genes.keys()
            newFamCount = len(panFams)
        else:
            newCoreFams = []
            for fam in fam2genes:
                if fam not in panFams:
                    panFams.append(fam)
                    newFamCount += 1
                if fam in coreFams:
                    newCoreFams.append(fam)
            coreFams = newCoreFams
        f.write("%s\t%d\t%d\t%d\t%d\n" %(sample, i+1, len(panFams), len(coreFams), newFamCount))
            
######## READ INPUT FILES ########
def readPslFile(file, best):
    #If best is True, only keep best matches. Otherwise keep all matches
    psls = {} #key = qName, val = Psl
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        psl = Psl(line)
        if psl.qName not in psls:
            psls[psl.qName] = [psl]
        elif best: #only take the best-matched transcripts
            currMatches = psls[psl.qName][-1].matches
            if psl.matches < currMatches:
                continue
            elif psl.matches == currMatches:
                psls[psl.qName].append(psl)
            else:
                psls[psl.qName] = [psl]
        else:
            psls[psl.qName].append(psl)
    f.close()
    return psls

def readSampleList(file):
    samples = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) > 0 and line[0] != '#':
            if line not in samples:
                samples.append(line)
    f.close()
    return samples

def readFasta(file):
    header2seq = {}
    f = open(file, 'r')
    header = ''
    seq = ''
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        if line[0] == '>': #done with old sequence, start a new one
            if header != '' and seq != '':
                if header in header2seq:
                    sys.stderr.write("Warning: repeatitive records of sequence %s in file %s\n" %(header, file))
                header2seq[header] = seq
            items = line.split()
            header = items[0].lstrip('>')
            seq = ''
        else:
            seq += line
   
    #last sequence
    if header != '' and seq != '':
        if header in header2seq:
            sys.stderr.write("Warning: repeatitive records of sequence %s in file %s\n" %(header, file))
        header2seq[header] = seq
    
    f.close()
    return header2seq

def readFasta2(file):
    header2seq = readFasta(file)
    fam2gene2seq = {}
    for header, seq in header2seq.iteritems():
        items = header.split(';')
        assert len(items) == 2
        fam = items[0]
        gene = items[1]
        if fam not in fam2gene2seq:
            fam2gene2seq[fam] = {gene:seq}
        else:
            fam2gene2seq[fam][gene] = seq
    return fam2gene2seq

def readFastas(files):
    header2seq = {}
    for file in files:
        h2s = readFasta(file)
        for h, s in h2s.iteritems():
            if h in header2seq:
                sys.stderr.write("Warning: repeatitive records of sequence %s in files %s\n" %(h,  ",".join(files)))
            header2seq[h] = s
    return header2seq

def writeFasta(header2seq, file):
    f = open(file, 'w')
    for header, seq in header2seq.iteritems():
        f.write(">%s\n%s\n" %(header, seq))
    f.close()

def writeFasta2(header2seq, fams, file):
    f = open(file, 'w')
    for fam in fams:
        #seq = header2seq[ fam.longestGene ]
        for g in fam.genes:
            seq = header2seq[ g ]
            header = "%s;%s" %(fam.id, g)
            f.write(">%s\n%s\n" %(header, seq))
    f.close()

def writeFasta3(fam2gene2seq, file):
    f = open(file, 'w')
    for fam, gene2seq in fam2gene2seq.iteritems():
        for gene, seq in gene2seq.iteritems():
            header = "%s;%s" %(fam, gene)
            f.write(">%s\n%s\n" %(header, seq))
    f.close()

####### MAIN ################
def main():
    usage = "%prog <fasta directory> <output file>"
    parser = OptionParser(usage = usage)
    parser.add_option('-s', '--sampleList', dest='sampleList', help='File containing list of samples to be included in the desired order. If not specified, includes all samples in the input directory')
    parser.add_option('-p', '--prot', dest='prot', action='store_true', default=False, help='Specified if input sequences are protein sequences. Default is for nucleotide sequences.')
    parser.add_option('-i', '--identity', dest='identity', type='float', default=90.0, help='Minimum alignment identity (%). Default=%default %')
    parser.add_option('-c', '--coverage', dest='coverage', type='float', default=90.0, help='Minimum alignment length coverage (%). Default=%default %')
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()

    if len(args) < 2:
        parser.error("Required 2 inputs\n")

    if options.sampleList:
        samples = readSampleList( options.sampleList )
    else:
        samples = os.listdir( args[0] )
    options.samples = samples

    i = Stack( Setup(args[0], samples, args[1], options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from referenceViz.src.pslUnion import *
    main()


