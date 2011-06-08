# libPlotting.py
# a library for repetitive plotting code in the assemblathon analysis
# dent earl, dearl (a) soe ucsc edu
#
# Library taken from assemblAnalysis developed by Dent Earl and others, 
# with modifications & additions by Ngan Nguyen to adapt to the reference project
#
##############################
# Copyright (C) 2009-2011 by 
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedict.paten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
# ... and other members of the Reconstruction Team of David Haussler's 
# lab (BME Dept. UCSC).
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
from numpy import *
from matplotlib.ticker import LogLocator

def initOptions( parser ):
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output, if --outFormat is all or png. default=%default')
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]. default=%default' )
   parser.add_option( '--out', dest='out', default='myPlot',
                      type='string',
                      help='filename where figure will be created. No extension needed. default=%default' )

def checkOptions( options, parser ):
   if options.dpi < 72:
      parser.error('I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.' % options.dpi )
   if options.outFormat not in ('pdf', 'png', 'eps', 'all'):
      parser.error('Unrecognized output format: %s. Choose one from: pdf png eps all.' % options.outFormat )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]

def initImage( width, height, options ):
   """
   initImage takes a width and height and returns
   both a fig and pdf object. options must contain outFormat,
   and dpi
   """
   import matplotlib.backends.backend_pdf as pltBack
   import matplotlib.pyplot as plt
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=( width, height ), dpi=options.dpi, facecolor='w' )
   return ( fig, pdf )

def writeImage( fig, pdf, options ):
   """
   writeImage assumes options contains outFormat and dpi.
   """
   if options.outFormat == 'pdf':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
   elif options.outFormat == 'png':
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'all':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
      fig.savefig( options.out + '.eps', format='eps' )
   elif options.outFormat == 'eps':
      fig.savefig( options.out + '.eps', format='eps' )


#Ngan's addition
def getColors( num ):             
    #Generate list of colors
    dimsize = (num + 2) ** (1/3.0)                                                            
    if int(dimsize) < dimsize:                                                          
        dimsize = int(dimsize) + 1                                                      
    else:                                                                               
        dimsize = int(dimsize)                                                          
    #dimsize += 2                                                                        
    colors = []                                                                         
                                                                                        
    vals = linspace( 0, 1, dimsize )                                                    
    for r in vals:                                                                      
        for g in vals:                                                                  
            for b in vals:                                                              
                colors.append( (r, g, b) )                                              
    return colors                                                                       
    
def getColors2( num ):
    #Generate list of colors and try to maximize the range these they cover
    colors = getColors( num )
    if len(colors) < num + 2:
        sys.stderr.write("colors array has less than required number of colors.\n") 
        sys.exit( 1 ) 
    indices = [1, len(colors) - 2]
    while  len(indices) < num:
        i = 0
        while i < len(indices) -1 and len(indices) < num:
            midpoint = (indices[i] + indices[i + 1])/2
            if midpoint != indices[i] and midpoint != indices[i+1]:
                indices.insert(i + 1, midpoint)
                i += 1
            i += 1
    colors2 = []
    for index in indices:
        colors2.append( colors[index] )
    return colors2


def setAxes( fig ):                                                                     
    return fig.add_axes( [0.12, 0.1, 0.83, 0.85] )                                      
                                                                                        
def editSpine( axes ):                                                                 
    for loc, spine in axes.spines.iteritems():                                          
        if loc in [ 'left', 'bottom' ]:                                                 
            spine.set_position( ('outward', 10) )                                       
        elif loc in ['right', 'top']:                                                   
            spine.set_color('none')                                                     
        else:                                                                           
            raise ValueError( 'Unknown spine location %s\n' % loc )                     
                                                                                        
def setTicks( axes ):                                                                   
    axes.xaxis.set_ticks_position( 'bottom' )                                           
    axes.yaxis.set_ticks_position( 'left' )                                             
    minorLocator = LogLocator( base=10, subs = range(1, 10) )                           
    axes.xaxis.set_minor_locator( minorLocator )                                        
                                                   


