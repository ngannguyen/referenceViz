#!/usr/bin/env python

#Wed Apr 10 15:30:53 PDT 2013
#
#Generates necessary files to make assembly hub
#Input: 1/ Output directory (e.g ~nknguyen/public_html/ecoli/hub/testhubs)
#       2/ hal file of the multiple alignment
#       3/ (list of samples)
#Output:
#   outdir/
#       hub.txt
#       genomes.txt
#

import os, sys, re
from optparse import OptionParser
from sonLib.bioio import system  #(may need to remove this dependency on sonLib)

def getLongestSeq(seq2len):
    seqs = sorted( [(seq, len) for seq, len in seq2len.iteritems()], key=lambda item:item[1], reverse=True )
    return seqs[0]

def getGenomeSequencesFromHal(halfile, genome):
    statsfile = "%s-seqStats.txt" %genome
    system("halStats --sequenceStats %s %s > %s" %(genome, halfile, statsfile))
    
    seq2len = {}
    f = open(statsfile, 'r')
    for line in f:
        if len(line) < 2 or re.search("SequenceName", line):
            continue
        items = line.strip().split(", ")
        seq = items[0].split('.')[-1]
        #seq = items[0]
        l = int(items[1])
        seq2len[seq] = l
    f.close()
    system("rm %s" %statsfile)

    return seq2len

def getGenomesFromHal(halfile):
    #Get a list of all genomes from the output of halStats
    statsfile = "halStats.txt"
    system("halStats --genomes %s > %s" %(halfile, statsfile))
    
    f = open(statsfile, 'r')
    genomes = f.readline().strip().split(',')
    f.close()

    #clean up
    system("rm %s" %statsfile)
    
    return genomes

def makeTwoBitSeqFile(genome, halfile, outdir):
    fafile = os.path.join(outdir, "%s.fa" %genome)
    system("hal2fasta --outFaPath %s %s %s" %(fafile, halfile, genome))
    
    #if sequence headers have "." (e.g genome.chr), reformat the header to only have "chr"
    fafile2 = "%s2" %fafile
    cmd = "awk '{ if($0 ~/>/){split($1, arr, \".\"); if(length(arr) > 1 ){print \">\" arr[2]}else{print $0} }else{ print $0} }' %s > %s" %(fafile, fafile2)
    system(cmd)
    system("rm %s" %fafile)

    #convert to 2bit files
    twobitfile = os.path.join(outdir, "%s.2bit" %genome)
    system("faToTwoBit %s %s" %(fafile2, twobitfile))
    system("rm %s" %fafile2)

def writeTrackDbFile(genomes, halfile, outdir):
    #Make softlink of halfile in outdir:
    #localHalfile = os.path.join(outdir, os.path.basename(halfile))
    #system("ln -s %s %s" %(os.path.abspath(halfile), localHalfile))
    
    filename = os.path.join(outdir, "trackDb.txt")
    f = open(filename, 'w')
    for genome in genomes:
        if re.search(genome, outdir): #current genome
            continue
        #SNAKE TRACKS
        f.write("track snake%s\n" %genome)
        f.write("longLabel %s Snake\n" %genome)
        f.write("shortLabel %s\n" %genome)
        f.write("otherSpecies %s\n" %genome)
        #f.write("visibility hide\n")
        f.write("visibility pack\n")
        f.write("bigDataUrl %s\n" % halfile)
        f.write("type halSnake\n")
        f.write("\n")
    f.close()

def writeDescriptionFile(genome, outdir):
    filename = os.path.join(outdir, "description.html")
    f = open(filename, 'w')
    f.write("%s\n" %genome)
    f.close()
    return 

def getLodFiles(halfile, options, outdir):
    lodtxtfile = os.path.join(outdir, "lod.txt") #outdir/lod.txt
    loddir = os.path.join(outdir, "lod") #outdir/lod
    if options.lodtxtfile and options.loddir: #if lod files were given, then just make soft links to them
        if os.path.exists(lodtxtfile):
            if os.path.abspath(lodtxtfile) != os.path.abspath(options.lodtxtfile):
                system("rm %s" %lodtxtfile)
                system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))

        if os.path.exists(loddir):
            if os.path.abspath(loddir) != os.path.abspath(options.loddir):
                if os.path.islink(loddir):
                    system("rm %s" %loddir)
                else:
                    system("rm -Rf %s" %loddir)
                loddir = os.path.join(outdir, os.path.basename(options.loddir))
                system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
    else: #if lod files were not given, create them using halLodInterpolate.py
        system("halLodInterpolate.py %s %s --outHalDir %s" %(halfile, lodtxtfile, loddir))
    return lodtxtfile, loddir

def writeGenomesFile(genomes, halfile, options, outdir):
    '''Write genome for all samples in hal file
    '''
    #Create lod files if useLod is specified
    lodtxtfile = ''
    loddir = ''
    if options.lod:
        lodtxtfile, loddir = getLodFiles(halfile, options, outdir)
    
    localHalfile = os.path.join(outdir, os.path.basename(halfile))
    system("ln -s %s %s" %(os.path.abspath(halfile), localHalfile))

    filename = os.path.join(outdir, "genomes.txt")
    f = open(filename, 'w')
    for genome in genomes:
        genomedir = os.path.join(outdir, genome)
        system("mkdir -p %s" % genomedir)
        f.write("genome %s\n" %genome)

        #create trackDb for the current genome:
        if lodtxtfile == '':
            writeTrackDbFile(genomes, "../%s" % os.path.basename(halfile), genomedir)
        else:
            writeTrackDbFile(genomes, "../%s" % os.path.basename(lodtxtfile), genomedir)
        f.write("trackDb %s/trackDb.txt\n" %genome)
        
        #create 2bit file for the current genome:
        makeTwoBitSeqFile(genome, halfile, genomedir)
        f.write("twoBitPath %s/%s.2bit\n" % (genome, genome))

        #other info
        f.write("groups groups.txt\n")

        writeDescriptionFile(genome, genomedir)
        f.write("htmlPath %s/description.html\n" %genome)
        f.write("organism %s\n" %genome)
        f.write("orderKey 4800\n")
        f.write("scientificName %s\n" %genome)
        
        seq2len = getGenomeSequencesFromHal(halfile, genome)
        (seq, l) = getLongestSeq(seq2len)
        f.write("defaultPos %s:1-%d\n" %(seq, l))
        f.write("\n")
    f.close()

def writeGroupFile(outdir):
    filename = os.path.join(outdir, "groups.txt")
    f = open(filename, 'w')
    f.write("name user\n")
    f.write("label Custom\n")
    f.write("priority 1\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")

    f.write("name map\n")
    f.write("label Mapping\n")
    f.write("priority 2\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")

    f.write("name x\n")
    f.write("label Experimental\n")
    f.write("priority 10\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")
    f.close()

def writeHubFile(outdir, options):
    hubfile = os.path.join(outdir, "hub.txt")
    f = open(hubfile, "w")
    f.write("hub %s\n" %options.hubLabel)
    f.write("shortLabel %s\n" %options.shortLabel)
    f.write("longLabel %s\n" %options.longLabel)
    f.write("genomesFile genomes.txt\n")
    f.write("email %s\n" %options.email)
    f.close()

def addOptions(parser):
    parser.add_option('--lod', dest='lod', action="store_true", default=False, help='If specified, create "level of detail" (lod) hal files and will put the lod.txt at the bigUrl instead of the original hal file. Default=%default')
    parser.add_option('--lodtxtfile', dest='lodtxtfile', help='"hal Level of detail" lod text file. If specified, will put this at the bigUrl instead of the hal file. Default=%default')
    parser.add_option('--loddir', dest='loddir', help='"hal Level of detail" lod dir. If specified, will put this at the bigUrl instead of the hal file. Default=%default')
    
    parser.add_option('--hub', dest='hubLabel', default='myHub', help='a single-word name of the directory containing the track hub files. Not displayed to hub users. Default=%default')
    parser.add_option('--shortLabel', dest='shortLabel', default='my hub', help='the short name for the track hub. Suggested maximum length is 17 characters. Displayed as the hub name on the Track Hubs page and the track group name on the browser tracks page. Default=%default')
    parser.add_option('--longLabel', dest='longLabel', default='my hub', help='a longer descriptive label for the track hub. Suggested maximum length is 80 characters. Displayed in the description field on the Track Hubs page. Default=%default')
    parser.add_option('--email', dest='email', default='NoEmail', help='the contact to whom questions regarding the track hub should be directed. Default=%default')
    #parser.add_option('')

def checkOptions(parser, args, options):
    if len(args) < 2:
        parser.error("Required two input arguments, %d was provided\n" %len(args))
    if not os.path.exists(args[0]):
        parser.error("Input hal file %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        system("mkdir -p %s" %args[1])
    elif not os.path.isdir(args[1]):
        parser.error("Output directory specified (%s) is not a directory\n" %args[1])

def main():
    usage = '%prog <halFile> <outputDirectory> [options]'
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    halfile = args[0]
    outdir = args[1]

    writeHubFile(outdir, options)
    writeGroupFile(outdir)
    genomes = getGenomesFromHal(halfile)
    
    writeGenomesFile(genomes, halfile, options, outdir)

if __name__ == '__main__':
    main()

