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

def getPdf( name ):
   import matplotlib.backends.backend_pdf as pltBack
   pdf = pltBack.PdfPages( name + '.pdf' )
   return pdf

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

def getColors0():
    colors = [ "#1f77b4", "#aec7e8", # blues 
               "#9467bd", "#c5b0d5", # lavenders
               "#ff7f0e", "#ffbb78", # oranges
               "#2ca02c", "#98df8a", # greens
               "#d62728", "#ff9896", # reds
               "#4A766E", "#48D1CC", #
               "#DAA520", "#DBDB70", #gold
               "#302B54", "#5959AB", #purple - blueish
               "#F6A8B6", "#F6CCDA", #pink
               "#7B68EE", "#6600FF",
               "#AA6600", "#B5A642",
               "#BC8F8F", "#C76114"] 

    return colors

def getColors1():
    colors = ["#E41A1C", #ref/hg19
              "#377EB8", "#4DAF4A", "#984EA3", #apd, cox, dbb
              "#FF7F00", "#FFFF33", "#A65628", #mann, mcf, qbl
              "#F781BF", #ssto
              "#C2E699", "#78C679", "#238443", #NA12878, NA12891, NA12892
              "#FED98E", "#FE9929", "#CC4C02", #NA19239, NA19238, NA19240
              "#7B68EE", "#6600FF", "#AA6600", "#B5A642", "#BC8F8F", "#C76114",
              "#999999"] #panTro
   
    return colors

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
                                                
def bihist( y1, y2, axes, bins, orientation, color=None ):
    n1, bins1, patch1 = axes.hist( y1, bins=bins, orientation=orientation, color=color ) #Top hist
    n2, bins2, patch2 = axes.hist( y2, bins=bins, orientation=orientation, color=color ) #bottom hist
    #set ymax:
    if orientation == 'vertical':
        ymax = max( [ i.get_height() for i in patch1 ] )
        for i in patch2:
            i.set_height( -i.get_height() )
        ymin = min( [ i.get_height() for i in patch2 ] )
    elif orientation == 'horizontal':
        ymax = max( [ i.get_width() for i in patch1 ] )
        for i in patch2:
            i.set_width( -i.get_width() )
        ymin = min( [ i.get_width() for i in patch2 ] )
    
    #axes.set_ylim( ymin*1.1, ymax*1.1 )
    return ymin, ymax

