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

def renameFaFile(filename, outdir, newseqname, old2new):
    currfilename = os.path.basename(filename)
    filenameItems = currfilename.split('.')
    
    if currfilename.isalnum():
        newfilename = currfilename
    elif len(filenameItems) > 1:
        newfilename = '.'.join(filenameItems[:-1])
        newfilename = "%s.%s" %(convertToAlnum(newfilename), filenameItems[-1])
    else:
        newfilename = convertToAlnum(currfilename)
    
    if currfilename != newfilename:
        old2new[currfilename] = newfilename
    newfilepath = os.path.join(outdir, newfilename)

    ifh = open(filename, 'r')
    ofh = open(newfilepath, 'w')
    count = 0
    for line in ifh:
        if len(line) > 1 and line[0] == '>':
            count += 1
            currheader = line.strip().lstrip('>')
            if newseqname:
                seqname = "%sSEQ%d" %(newfilename.split('.')[0], count)
            else:
                items = currheader.split()
                subitems = items[0].rstrip('|').split('|')
                seqname = convertToAlnum(subitems[-1])
            ofh.write(">%s\n" %seqname)

            if seqname != currheader:
                old2new[currheader] = seqname
        else:
            ofh.write(line)
    ifh.close()
    ofh.close()
    return

def renameInputs(input, outdir, newseqname):
    old2new = {}
    if os.path.isdir(input): # directory
        #Rename directory
        currdirname = os.path.basename(input)
        newdirname = currdirname
        if not currdirname.isalnum():
            newdirname = convertToAlnum( currdirname )
            old2new[currdirname] = newdirname

        newdir = os.path.join( outdir, newdirname )
        system( "mkdir -p %s" % newdir )
        
        #Going through the files in the directory and rename them and their sequences
        files = os.listdir(input)
        for file in files:
            renameFaFile( os.path.join(input, file), newdir, newseqname, old2new )
    else: #file
        renameFaFile( input, outdir, newseqname, old2new )
    
    #Print file that maps old names to new names:
    mapfile = os.path.join(outdir, 'headers_old2new.txt')
    if os.path.exists(mapfile):
        system('rm %s' %mapfile)
    f = open(mapfile, 'w')
    for o in sorted(old2new.keys()):
        f.write("%s\t%s\n" %(o, old2new[o]))
    f.close()
    return

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
    parser.add_option('--new', dest='new', action='store_true', default=False, help='If specified, will rename the sequences of each file with <filenameSeq#> with # is the order the sequence appear in the file. Default=%default')
    options, args = parser.parse_args()
    checkInputs(parser, args, options)

    renameInputs(args[0], args[1], options.new)

if __name__ == '__main__':
    main()


