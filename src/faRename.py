#!/usr/bin/env python

"""
nknguyen at soe ucsc edu
Fri Jan 25 16:09:44 PST 2013

Input: a fasta file or a directory of fasta files
       output directory
Rename the directories, the file names, and the sequence names so that they only contain alphanumeric values
"""
import os, sys, re, time
from sonLib.bioio import system
from optparse import OptionParser

def convertToAlnum(inputStr):
    items = re.split("[^a-zA-Z0-9]", inputStr)
    items = [item.capitalize() for item in items] 
    newStr = "".join(items)
    assert newStr.isalnum()
    return newStr

def renameFaFile(filename, outdir):
    filenameItems = os.path.basename(filename).split('.')
    newfilename = '.'.join(filenameItems[:-1])
    if len(filenameItems) > 1:
        newfilename = os.path.join(outdir, "%s.%s" %(convertToAlnum(newfilename), filenameItems[-1]))
    else:
        newfilename = os.path.join(outdir, convertToAlnum(newfilename))

    ifh = open(filename, 'r')
    ofh = open(newfilename, 'w')
    for line in ifh:
        if len(line) > 1 and line[0] == '>':
            items = line.lstrip('>').split()
            subitems = items[0].rstrip('|').split('|')
            seqname = convertToAlnum(subitems[-1])
            ofh.write(">%s\n" %seqname)
        else:
            ofh.write(line)
    ifh.close()
    ofh.close()
    return

def renameInputs(input, outdir):
    if os.path.isdir(input): # directory
        #Rename directory
        dirname = convertToAlnum( os.path.basename(input) )
        newdir = os.path.join( outdir, dirname )
        system( "mkdir -p %s" % newdir )
        
        #Going through the files in the directory and rename them and their sequences
        files = os.listdir(input)
        for file in files:
            renameFaFile( os.path.join(input, file), newdir )
    else: #file
        renameFaFile(input, outdir)

def checkInputs(parser, args, options):
    if len(args) < 2:
        parser.error("Please specify inputFile/Directory and outputDirectory\n")
    input = args[0]
    outdir = args[1]
    if not os.path.exists(input):
        parser.error("Input %s does not exist\n" %input)
    if not os.path.exists(outdir):
        system("mkdir %s" %outdir)
    return

def main():
    usage = ('Usage: %prog inputFile/Directory outputDirectory \n')
    parser = OptionParser( usage = usage )
    options, args = parser.parse_args()
    checkInputs(parser, args, options)

    renameInputs(args[0], args[1])

if __name__ == '__main__':
    main()


