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
    f.write("\\begin{tabular}{c|c|r|r|r}\n")
    f.write("\\multicolumn{5}{c}{%s} \\\\\n" %title)
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("\\multicolumn{2}{c|}{Sample} & Deletion & Non-linear Operation & SNP\\\\\n")
    #f.write("\\multicolumn{2}{c}{\\multirow{2}{*}{}} & \\multicolumn{2}{c}{Conserved} & \\multicolumn{2}{c}{Non-conserved} & \\multicolumn{2}{c}{Total} \\\\\n")
    #f.write("\\cline{3-8}\n")
    #f.write("\\multicolumn{2}{c}{} & Count & Percentage & Count & Percentage & Count & Percentage\\\\\n")
    f.write("\\hline\n")

def prettyFloat( num ):
    if num == 0:
        return "0"
    else:
        return "%.2e" %num

def tab( f, samplesList, sampleNames ):
    refname1 = samplesList[0][0].attrib['referenceName']
    refname2 = samplesList[1][0].attrib['referenceName']
    for s in sampleNames:
        #altColor = 1
        for altColor in [1,0]:
            #Get #Deletions, #Non-linearOps
            numDels = -1
            numDelsPerAlignedBase = -1
            numNonLinearOps = -1
            numNonLinearOpsPerAlignedBase = -1
            for sample in samplesList[altColor]:
                if sample.attrib['sampleName'] == s:
                    numDels = int( sample.attrib['totalDeletion'] )
                    numDelsPerAlignedBase = float( sample.attrib['totalDeletionPerAlignedBase'] )
                    numDelsPerAlignedBase = prettyFloat(numDelsPerAlignedBase)

                    #numNonLinearOps = int(sample.attrib['totalIntraJoin']) + int(sample.attrib['totalInterJoin'])
                    #numNonLinearOpsPerAlignedBase = float(sample.attrib['totalInterJoinPerAlignedBase']) + float(sample.attrib['totalIntraJoinPerAlignedBase'])
                    numNonLinearOps = int(sample.attrib['totalIntraJoin'])
                    numNonLinearOpsPerAlignedBase = float(sample.attrib['totalIntraJoinPerAlignedBase'])
                    numNonLinearOpsPerAlignedBase = prettyFloat(numNonLinearOpsPerAlignedBase)
                    break
            #Get the Snps#
            numSnps = -1
            numSnpsPerAlignedBase = -1
            for sample in samplesList[altColor + 2]:
                if sample.attrib['sampleName'] == s:
                    numSnps = int(sample.attrib['totalErrors'])
                    numSnpsPerAlignedBase = numSnps/float( sample.attrib['totalCalls'])
                    numSnpsPerAlignedBase = prettyFloat( numSnpsPerAlignedBase )
                    break
           
            if altColor == 1:
                f.write("\\multirow{2}{*}{%s} &\\cellcolor[gray]{0.9} %s & \\cellcolor[gray]{0.9} %d (%s) & \\cellcolor[gray]{0.9} %d (%s) & \\cellcolor[gray]{0.9} %d (%s) \\\\\n" % \
                        (s, refname2, numDels, numDelsPerAlignedBase, numNonLinearOps, numNonLinearOpsPerAlignedBase, numSnps, numSnpsPerAlignedBase))
            else:
                f.write("& %s & %d (%s) & %d (%s) & %d (%s) \\\\\n" %\
                        (refname1, numDels, numDelsPerAlignedBase, numNonLinearOps, numNonLinearOpsPerAlignedBase, numSnps, numSnpsPerAlignedBase))
                f.write("\\hline\n\n")
            #altColor = 1 - altColor
    #f.write("\\hline\n\n")

def tableCloser(f, captionStr, label):
    f.write("\\end{tabular}\n")
    f.write("}\n")
    f.write("\\caption{%s}\n" %captionStr)
    f.write("\\end{center}\n")
    f.write("\\label{%s}" %label)
    f.write("\\end{table}\n\n")

def makeLatexTab( samplesList, outdir, sampleNames ):
    #samples = sorted( samples, key=lambda s: s.attrib[ 'sampleName' ] )
    absOutdir = os.path.abspath( outdir )
    if not os.path.exists( absOutdir ):
        sys.stderr.write( 'Output directory does not exist\n' )
        sys.exit( 1 )
    if len(sampleNames) < 1:
        return

    refname1 = samplesList[0][0].attrib['referenceName']
    refname2 = samplesList[1][0].attrib['referenceName']
    
    if refname1 != samplesList[2][0].attrib['referenceName'] or refname2 != samplesList[3][0].attrib['referenceName']:
        sys.stderr.write("Input file 1 and input file 3 must have the same referenceName;\
                          and input file 2 and input file 4 must have the same referenceName.")
        sys.exit( 1 )
    
    out = os.path.join( absOutdir, "indelStats_%s-%s.tex" % (refname1,refname2) )
    f = open( out, 'w' )

    writeDocumentStart( f )
    title = "`Errors' Statistics"
    tabHeader( f, title )
    tab( f, samplesList, sampleNames )
    captionStr = "%s-%s"    % (refname1,refname2)
    label = "%s-%s" %(refname1, refname2)
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
    usage = "Usage %prog [options] pathStats*.xml1 pathStats*.xml2 snpStats*.xml1 snpStats*.xml2"
    parser = OptionParser( usage = usage )
    parser.add_option('--outdir', dest='outdir', default='.', help="Output directory")
    parser.add_option('--samples', dest='samples', default='apd,cox,dbb,mann,mcf,qbl,ssto,NA12891,NA12892,NA12878,NA19239,NA19238,NA19240,panTro2', help="Comma separated list of samples in the desired order")
    options, args = parser.parse_args()
    
    if len(args) < 4:
        parser.error("Please specify two pathStatsXml files and two snpStats Files\n")
    samplesList = readfiles( args ) #xml trees
    sampleNames = options.samples.split(',')
    makeLatexTab( samplesList, options.outdir, sampleNames )


if __name__ == '__main__':
    main()
