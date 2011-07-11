#!/usr/bin/env python

"""
nknguyen at soe ucsc edu
Jun 29 2011
This script generates summary statistics for different input parameters of (cactus) reference algorithm
"""

import os, sys, re, time
import xml.etree.ElementTree as ET
from optparse import OptionParser

#from numpy import *
from numpy import array
import libPlotting as libplot
import matplotlib.pyplot as pyplot
from matplotlib.ticker import *
from matplotlib.font_manager import FontProperties

class Run():
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename( path )
        self.params = None
        self.sortOrder = None
        self.data = {}

    def setParams(self, patternStr, p2type):
        m = re.search( patternStr, self.name )
        if m == None:
            sys.stderr.write( "Directory '%s' does not match pattern '%s'\n" %(self.name, patternStr) )
            sys.exit( 1 )
        self.params = m.groupdict()
        self.p2type = p2type

    def setSortOrder(self, sortOrder):
        self.sortOrder = sortOrder.split(',') 

    def setData(self, analysisName, refname, data, colNames):
        if analysisName not in self.data:
            self.data[ analysisName ] = {}

        self.data[ analysisName ][ refname ] = [data, colNames]

    def printRun(self, filehandle):
        filehandle.write("Path: %s\n" % self.path)
        filehandle.write("Name: %s\n" % self.name)
        if self.params:
            filehandle.write("Params:\t" )
            for p in self.params:
                filehandle.write( "%s:%s\t" %(p, self.params[p]))
            filehandle.write("\n")
        if self.sortOrder:
            filehandle.write("SortOrder: %s \n" %(','.join(self.sortOrder) ))

    def __cmp__(self, otherRun):
        params = self.sortOrder
        if len(params) == 0:
            return cmp(self.name, otherRun.name)

        i = 0
        while i < len(params):
            currParam = params[i]
            
            v = self.params[ currParam ]
            vother = otherRun.params[ currParam ]
            if self.p2type[ currParam ] == 'integer':
                v = int( v )
                vother = int( vother )
                #print "Found integers: %d, %d" % (v, vother)
                #print cmp(v, vother)

            currCmp = cmp( v, vother )
            if currCmp == 0:
                i += 1
            else:
                return currCmp
        return 0

def getRuns( indir, patternStr, params ):
    runs = []
    p2type = {}
    for p in params:
        p2type[p.attrib['name']] = p.attrib['type']
    #print p2type

    for file in os.listdir( indir ):
        if not os.path.isdir( os.path.join(indir, file) ):
            continue
        run = Run( os.path.join(indir, file) )
        run.setParams( patternStr, p2type )
        runs.append( run )
    return runs            

def readConfig( file ):
    xmltree = ET.parse( file )
    root = xmltree.getroot()
    params = root.find( 'parameters' )
    if params == None:
        sys.stderr.write( 'No parameters were specified in the config file!\n' )
        sys.exit( 1 )

    params = params.findall( 'param' )

    #Sorting orders:
    sortOrders = root.find( 'sortOrders' )
    if sortOrders:
        sortOrders = [ so.attrib[ 'order' ] for so in sortOrders.findall('sortOrder') ]
    
    #if no sortOrder was specified in the config file, use the params order
    if sortOrders == None or len(sortOrders) == 0: 
        names = [ param.attrib[ 'name' ] for param in params ]
        sortOrders = [ ','.join(names) ]
    
    #Analyses to look at:
    analyses = root.find( 'analyses' )
    if analyses:
        analyses = analyses.findall( 'analysis' )
    if analyses == None or len(analyses) == 0:
        sys.stderr.write( 'No analyses was specified in the config file. Nothing to do, exit!\n' )
        sys.exit( 0 )

    return (params, sortOrders, analyses)

def getDirNamePattern( params ):
    pattern = []
    for param in params:
        vals = param.attrib[ 'values' ].replace(',', '|')
        pattern.append( '(?P<%s>%s)' %(param.attrib[ 'name' ], vals) )
    patternStr =  '-'.join( pattern ) + '$'
    return patternStr

def getData2( runs, analysisName, refname ):
    rows = []
    rowNames = []
    colNames = []
    for run in runs:
        if analysisName not in run.data:
            continue
        #print "%s, %s, %s" % (run.name, analysisName, refname)
        rowNames.append( run.name )
        #print run.data
        rundata = run.data[ analysisName ]
        #print rundata
        items = rundata[ refname ]
        #print items
        rows.append(items[0])
        if len(colNames) == 0:
            colNames = items[1]
    
    return array( rows ), colNames, rowNames

def getData( indir, runs, infile, sortField, keyField, analysisName ):
    rows = []
    sampleNames = []
    refname = ''
    #print "getData..."
    for run in runs:
        row = []
        file = os.path.join(indir, run.name, infile)
        if not os.path.exists( file ):
            continue
        #print file
        xmltree = ET.parse( file )
        root = xmltree.getroot()
        allsamples = root.findall('statsForSample')
        samples = []
        for s in allsamples:
            if s.attrib[ 'sampleName' ] != '' and s.attrib[ 'sampleName' ] != 'ROOT' and (s.attrib['sampleName'] not in ['hg19', 'reference']):
                samples.append( s )
        samples = sorted( samples, key=lambda sample:sample.attrib[sortField] )
        
        if len(sampleNames) == 0:
            refname = samples[0].attrib['referenceName']
            for s in samples:
                sampleNames.append( s.attrib['sampleName'] )
            sampleNames.append( 'Average' )
            #print "Refname: %s" %refname
            #print "sampleNames %s" %sampleNames
        
        for s in samples:
            row.append( float(s.attrib[keyField])  )
            #row.append( float(s.attrib['correctPerAligned'])  )
        avg = sum(row)/len(row)
        row.append(avg)
        run.setData(analysisName, refname, row, sampleNames)
        
    return refname

def setAxes( fig ):
    axesList = []
    axleft = 0.3
    axright = 0.98
    axwidth = axright - axleft
    axbottom = 0.08
    axtop = 0.95
    axheight = axtop - axbottom
    margin = 0.01

    plotw = (axwidth - margin)/2.0
    axesList.append( fig.add_axes( [axleft, axbottom, plotw, axheight] ) ) #left plot
    axesList.append( fig.add_axes( [axleft + plotw + margin, axbottom, plotw, axheight] ) ) #right plot

    return axesList


def drawSum( fig, runs, analysis, options, refnames ):
    analysisName = analysis.attrib[ 'name' ]

    axesList = []
    if len(refnames) == 1:
        axesList.append( libplot.setAxes(fig) )
    else:
        axesList = setAxes( fig )

    sortField = 'sampleName'
    #fontP = FontProperties()
    #fontP.set_size('xx-small')
    for i in range( len(refnames) ):
        #data, colLabels, rowLabels, refname = getData( options.indir, runs, files[i], sortField )
        #dataFile = os.path.join( options.outdir, 'sumData_' + files[i].strip('.xml') + '_' + sortId)
        #if not os.path.exists( dataFile ):
        #    data, colLabels, rowLabels, refname = getData( options.indir, runs, files[i], sortField )
        #    writeData( dataFile, data, colLabels, rowLabels, refname )
        #else:
        #    data, colLabels, rowLabels, refname = readData( dataFile )
        
        data, colLabels, rowLabels = getData2( runs, analysisName, refnames[i] )
        axes = axesList[i]
        caxes = axes.imshow( data, interpolation= 'nearest' )
        axes.set_title( refnames[i] )
        #axes.set_xticks( range(0, len(colLabels)), colLabels, rotation=45, fontproperties=fontP )
        axes.set_xticks( range(0, len(colLabels)) )
        axes.set_xticklabels(colLabels)
        for label in axes.xaxis.get_ticklabels():
            label.set_rotation( 90 )
            label.set_fontsize(4)
            #label.fontproperties = fontP
        axes.set_yticks( range(len(runs)) )
        if i > 0:
            rowLabels = [ '' for l in rowLabels ]
        axes.set_yticklabels( rowLabels )
        for label in axes.yaxis.get_ticklabels():
            label.set_fontsize(4)
        #else:
        #    axes.yaxis.set_major_locator( NullLocator() )
        #    axes.yaxis.set_major_formatter( NullFormatter() )
    
    #cbar = fig.colorbar( caxes, ticks=[
    #cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
    #cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])# horizontal colorbar
        


def getSummary( runs, analysis, options, sortId, refnames ):
    analysisName = analysis.attrib[ 'name' ]
    options.out = os.path.join( options.outdir, 'summary_%s_%d' %(analysisName, sortId) )
    #fig, pdf = libplot.initImage( 8.0, 10.0, options )
    fig, pdf = libplot.initImage( 16.0, 20.0, options )
    
    drawSum( fig, runs, analysis, options, refnames)
    libplot.writeImage( fig, pdf, options )
    

def initOptions( parser ):
    parser.add_option('-i', '--indir', dest='indir', help='Required argument. Location of the input directory. This input directory should contain subdirectories, each of which corresponds to one reference run')
    parser.add_option('-o', '--output', dest='outdir', help='Required argument. Output directory.')
    parser.add_option('-c', '--config', dest='config', help='Required argument. Config file specifying accepted values of different parameters, as well as what order to sort the output stats')
    parser.add_option('-s', '--silent', dest='silent', action="store_true", help="If specified, will not print out plots." )

def checkOptions( options, parser ):
    if not options.indir:
        parser.error( 'Location of input direcotry is required and not given.\n' )
    if not os.path.isdir(options.indir):
        parser.error( "%s is not a directory \n" % options.indir )
    
    if not options.outdir:
        parser.error( 'Location of output directory is required and not given.\n' )
    if not options.config:
        parser.error( 'Config file is required and not given.\n' )
    
    params, sortOrders, analyses = readConfig( options.config )
    options.params = params
    options.sortOrders = sortOrders
    options.analyses = analyses
    return
    
def dumpData( runs, analysisName, outdir, refnames ):
    outfile = os.path.join( outdir, "_".join(["summary", analysisName]) + ".txt" )
    f = open( outfile, 'w' )

    params = (runs[0].params).keys()
    colNames = [ p for p in params ]

    i = 0
    for run in runs:
        if analysisName not in run.data:
            continue
        i += 1
        rowParams = []
        for param in params:
            rowParams.append( run.params[param] )
        
        rundata = run.data[ analysisName ]
        rowData = []
        for ref in refnames:
            items = rundata[ ref ]

            #if len(colNames) == len( params ):
            if i == 1:
                colNames.extend( items[1] )
        
            rowData.extend([ str(d) for d in items[0] ])
        if i == 1:
            f.write( "%s\n" %( '\t'.join(colNames) ) )

        f.write( "%s\t%s\n" %( '\t'.join(rowParams), '\t'.join(rowData) ) )
    f.close()

def main():
    usage = ( 'usage: %prog [options]'
              '%prog generates summary stats of reference runs of different parameters' )
    parser = OptionParser( usage = usage )
    initOptions( parser )
    libplot.initOptions( parser )
    
    options, args = parser.parse_args()
    checkOptions( options, parser )
    libplot.checkOptions( options, parser )
    
    dirPatternStr = getDirNamePattern( options.params )
    runs = getRuns( options.indir, dirPatternStr, options.params )

    for analysis in options.analyses:
        #Read data:
        sortField = 'sampleName'
        files = analysis.attrib[ 'files' ].split(',')
        analysisName = analysis.attrib[ 'name' ]
        if len(files) == 0:
            sys.stderr.write('No input file(s) for analysis %s\n' %(analysis.attrib['name']) )
            continue
        
        refnames = []
        for file in files:
            refname = getData( options.indir, runs, file, sortField, analysis.attrib['keyField'], analysis.attrib[ 'name' ] )
            refnames.append(refname)
       
        
        dumpData( runs, analysis.attrib[ 'name' ], options.outdir, refnames )
        if options.silent:
            return

        s = 1
        for sortOrder in options.sortOrders:
            #Sort run with 'sortOrder'
            for run in runs:
                run.setSortOrder( sortOrder )
            runs.sort()

            getSummary( runs, analysis, options, s, refnames )
            s += 1


if __name__ == '__main__':
    main()

