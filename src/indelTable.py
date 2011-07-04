#!/usr/bin/env python

import os, sys, re
from optparse import OptionParser
import xml.etree.ElementTree as ET


def writeDocumentStart( f ):
    f.write( "\\documentclass[11pt]{article}\n" )

    f.write( "\\usepackage{epsfig}\n" )
    f.write( "\\usepackage{multirow}\n" )
    f.write( "\\usepackage{graphicx}\n" )
    f.write( "\\usepackage{array}\n" )
    f.write( "\\usepackage{color}\n" )
    f.write( "\\usepackage{rotating}\n" )
    f.write( "\\usepackage[table]{xcolor}\n" )
    f.write( "\n" )

    f.write( "\\newcommand{\\figref}[1]{Figure~\\ref{fig:#1}}\n" )
    f.write( "\\newcommand{\\tabref}[1]{Table~\\ref{tab:#1}}\n" )
    f.write( "\n" )

    f.write("\\textwidth=6.5in\n")
    f.write("\\textheight=9in\n")
    f.write("\\oddsidemargin=0in\n")
    f.write("\\evensidemargin=0in\n")
    f.write("\\topmargin=0in\n")
    f.write("\\topskip=0in\n")
    f.write("\\headheight=0in\n")
    f.write("\\headsep=0in\n")
    f.write("\n")

    f.write("\\begin{document}\n")
    return

def writeDocumentEnd( f ):
    f.write( "\\end{document}\n" )

def tabHeader( f, title ):
    f.write("\\begin{table}\n")
    f.write("\\begin{center}\n")
    f.write("\\scalebox{1}{%\n")
    #f.write("\\begin{tabular}{|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|>{\small}c|}\n")
    f.write("\\begin{tabular}{ccccccc}\n")
    f.write("\\multicolumn{7}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("Sample & SequenceLength & Insertion & Deletion & Indel & IntraJoin & InterJoin\\\\\n")
    #f.write("\\multicolumn{2}{c}{\\multirow{2}{*}{}} & \\multicolumn{2}{c}{Conserved} & \\multicolumn{2}{c}{Non-conserved} & \\multicolumn{2}{c}{Total} \\\\\n")
    #f.write("\\cline{3-8}\n")
    #f.write("\\multicolumn{2}{c}{} & Count & Percentage & Count & Percentage & Count & Percentage\\\\\n")
    f.write("\\hline\n")

def tab( f, samples ):
    altColor = 1
    for s in samples:
        if altColor == 1:
            f.write("\\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %s \\\\\n" % \
                    (s.attrib['sampleName'], s.attrib['totalSampleGenomeLength'], s.attrib['totalInsertion'], s.attrib['totalDeletion'], s.attrib['totalInsertionAndDeletion'], s.attrib['totalIntraJoin'], s.attrib['totalInterJoin']))
        else:
            f.write("%s & %s & %s & %s & %s & %s & %s \\\\\n" %\
                    (s.attrib['sampleName'], s.attrib['totalSampleGenomeLength'], s.attrib['totalInsertion'], s.attrib['totalDeletion'], s.attrib['totalInsertionAndDeletion'], s.attrib['totalIntraJoin'], s.attrib['totalInterJoin']))
        altColor = 1 - altColor
    f.write("\\hline\n\n")

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def makeLatexTab( samples, outdir ):
    samples = sorted( samples, key=lambda s: s.attrib[ 'sampleName' ] )
    absOutdir = os.path.abspath( outdir )
    if not os.path.exists( absOutdir ):
        #system("mkdir -p %s" %absOutdir)
        sys.stderr.write( 'Output directory does not exist\n' )
        sys.exit( 1 )
    if len(samples) < 1:
        return
    refname = samples[0].attrib['referenceName']
    out = os.path.join( absOutdir, "indelStats_%s.tex" % refname )
    f = open( out, 'w' )

    writeDocumentStart( f )
    title = "`Errors' Statistics, %s" %refname
    tabHeader( f, title )
    tab( f, samples )
    captionStr = "%s"    % samples[0].attrib['referenceName']
    label = samples[0].attrib['referenceName']
    tableCloser(f, captionStr, label)
 
    writeDocumentEnd( f )
    f.close()

def readfile( file ):
    if not os.path.exists( file ):
        sys.stderr.write( 'File %s does not exist\n' % file )
        exit( 1 )
    return ET.parse( file )
    
def readfiles( files ):
    samplesList = []
    for f in files:
        tree = readfile( f )
        root = tree.getroot()
        samples = []
        for sample in root.findall( 'statsForSample' ):
            name = sample.attrib[ 'sampleName' ]
            if name != "ROOT" and name != "":
                samples.append( sample )
        samplesList.append( samples )

    return samplesList

def main():
    usage = "Usage %prog [options] pathStats*.xml1 [pathStats*.xml2 ...]"
    parser = OptionParser( usage = usage )
    parser.add_option('--outdir', dest='outdir', default='.', help="Output directory")
    options, args = parser.parse_args()

    samplesList = readfiles( args ) #xml trees
    for samples in samplesList:
        makeLatexTab( samples, options.outdir )


if __name__ == '__main__':
    main()
