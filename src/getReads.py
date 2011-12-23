#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import getTempFile

from sonLib.bioio import logger
from sonLib.bioio import system
from sonLib.bioio import nameValue
from sonLib.bioio import getTempDirectory
from sonLib.bioio import setLogLevel

#Download, Unzip files & organize them to the output directory in this format:
#Outdir/
#   Sample1/
#       ILLUMINA/
#           PAIRED/
#               file1_reads1
#               file1_reads2
#               file2_reads1
#               file2_reads2
#               ...
#           SINGLE/
#               file1
#               file2
#               ...
#       LS454/
#           PAIRED/
#           SINGLE/
class Setup(Target):
    def __init__(self, options):
        Target.__init__(self, time=0.00025)
        self.options = options

    def run(self):
        setLogLevel("INFO")
        file2info = readSeqIndex(self.options.seqIndexFile, self.options.samples)
        #logger.info("Done reading sequence.index file\n")
        for sample in file2info:
            self.addChildTarget( RunSample(self.options, file2info[sample]) )

class RunSample(Target):
    def __init__(self, options, file2info):
        Target.__init__(self, time=0.00005)
        self.options = options
        self.file2info = file2info

    def run(self):
        filename2info = self.file2info 
        for file in filename2info:
            readfile = filename2info[file]

            samplePath = os.path.join( self.options.outdir, readfile.sample, readfile.platform, readfile.libLayout )
            system("mkdir -p %s" %(samplePath))
            samplefile = os.path.join(samplePath, readfile.name)

            sampleAddress = "%s/%s" %(self.options.ftpaddress, readfile.path)
            if not os.path.exists( samplefile.rstrip('.gz') ):
                system( "ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty -Tr -Q -l 500M -L- %s %s" %(sampleAddress, samplePath) )

        #Done downloading all the files, now gunzip them:
        self.setFollowOnTarget( Uncompress(self.options, filename2info) )
        #self.addChildTarget( GetFile(self.options, readfile['sample'], readfile['platform'], readfile['libLayout'], readfile['name'], readfile['path']) )

class Uncompress(Target):     
    def __init__(self, options, file2info):
        Target.__init__(self, time=0.000001)
        self.options = options
        self.file2info = file2info

    def run(self):
        file2info = self.file2info
        for file in self.file2info:
            readfile = file2info[file]
            samplePath = os.path.join( self.options.outdir, readfile.sample, readfile.platform, readfile.libLayout )
            samplefile = os.path.join(samplePath, readfile.name)
            if not os.path.exists( samplefile.rstrip('.gz') ):
                self.addChildTarget( UncompressFile(samplefile) )

class UncompressFile(Target):     
    def __init__(self, samplefile):
        Target.__init__(self, time=0.1)
        self.file = samplefile

    def run(self):
        system("gunzip %s" %self.file)

#class GetFile(Target):
#    def __init__(self, options, sample, platform, libLayout, name, path):
#        Target.__init__(self, time=0.00025)
#        #Target.__init__(self)
#        self.options = options
#        self.sample = sample
#        self.platform = platform
#        self.libLayout = libLayout
#        self.name = name
#        self.path = path
#
#    def run(self):
#        #logger.info("GETFILE\n")
#        samplePath = os.path.join( self.options.outdir, self.sample, self.platform, self.libLayout )
#        system("mkdir -p %s" %(samplePath))
#        samplefile = os.path.join(samplePath, self.name)
#
#        sampleAddress = "%s/%s" %(self.options.ftpaddress, self.path)
#        if not os.path.exists( samplefile.rstrip('.gz') ):
#            system( "ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty -Tr -Q -l 500M -L- %s %s" %(sampleAddress, samplePath) )
#            #system("wget -O - %s > %s" %(sampleAddress, samplefile))
#            system("gunzip %s" %samplefile)

class Readfile():
    def __init__(self, line):
        items = line.strip().split('\t')
        l = items[0].split('/')
        self.path = items[0]
        self.name = l[ len(l) -1 ]
        self.sample = items[9]
        self.platform = items[12].lower() #ILLUMINA, LS454, ABI_SOLID
        self.insertSize = items[17] 
        self.libLayout = items[18].lower() #PAIRED, OR SINGLE
        self.pairedFile = items[19]
        self.withdrawn = items[20]
        self.withdrawnDate = items[21]

        #if self.pairedFile == "":#If paired, but not mate file, then set as single
        #    self.libLayout = "single"
        self.readCount = items[23]
        self.baseCount = items[24]

#def readLine(line):
#    readfile = {}
#    items = line.strip().split('\t')
#    l = items[0].split('/')
#    readfile['path'] = items[0]
#    readfile['name'] = l[ len(l) -1 ]
#    readfile['sample'] = items[9]
#    readfile['platform'] = items[12].lower() #ILLUMINA, LS454, ABI_SOLID
#    readfile['insertSize'] = items[17] 
#    readfile['libLayout'] = items[18].lower() #PAIRED, OR SINGLE
#    readfile['pairedFile'] = items[19]
#    readfile['withdrawn'] = items[20]
#    readfile['withdrawnDate'] = items[21]
#
#    #if self.pairedFile == "":#If paired, but not mate file, then set as single
#    #    self.libLayout = "single"
#    readfile['readCount'] = items[23]
#    readfile['baseCount'] = items[24]
#    return readfile

def readSeqIndex(file, samples):
    #Read sequence.index file:
    infoFh = open(file, 'r')
    file2info = {} #key = sample, val = {filename:Readfile}

    for line in infoFh.readlines():
        rf = Readfile( line )
        if rf.sample not in samples:
            continue
        elif rf.pairedFile == "" and rf.libLayout == "paired":
            continue
        elif rf.withdrawn == '1' or rf.withdrawnDate != '':
            continue
        elif rf.sample not in file2info:
            file2info[ rf.sample ] = { rf.name: rf }
        else:
            file2info[ rf.sample ][ rf.name ] = rf
        #if rf['sample'] not in samples:
        #    continue
        #elif rf['pairedFile'] == "" and rf['libLayout'] == "paired":
        #    continue
        #elif rf['withdrawn'] == '1' or rf['withdrawnDate'] != '':
        #    continue
        #elif rf['sample'] not in file2info:
        #    file2info[ rf['sample'] ] = { rf['name'] : rf }
        #elif rf['name'] not in file2info[ rf['sample'] ]:
        #    file2info[ rf['sample'] ][rf['name']] = rf
    infoFh.close()
    return file2info

def checkOptions(options, args, parser ):
    if len(args) < 4:
        parser.error("Need 3 input files, only %d provided\n" %(len(args)))
    if not os.path.exists(args[0]):
        parser.error("sequence.index file %s does not exist\n" %(args[0]))
    system("mkdir -p %s" %args[2])
    options.seqIndexFile = args[0]
    options.samples = args[1].split(',')
    options.outdir = args[2]
    options.ftpaddress = args[3] #fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp

def main():
    #ftpAddress: fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp
    usage = "Usage: %prog [options] <sequence.index> <sample1[,sample2,...]> <outdir> <ftpAddress>"
    parser = OptionParser( usage = usage )
    Stack.addJobTreeOptions(parser)

    #initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( options, args, parser )

    i = Stack( Setup(options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobTree contains %d failed jobs\n" %i)

if __name__ == "__main__":
    from referenceViz.src.getReads import *
    main()

