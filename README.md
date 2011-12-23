#Reference MHC Project Analysis Codebase

This is a set of scripts developed by [Ngan Nguyen ](https://github.com/ngannguyen/) and [Benedict Paten ](https://github.com/benedictpaten/) to generate Tables and Figures for the Reference MHC project. 

##Dependencies
* [jobTree](https://github.com/benedictpaten/jobTree)
* [matplotlib ](http://matplotlib.sourceforge.net/)
* [sonLib](https://github.com/benedictpaten/sonLib)
* python version 2.6 or after, but before 3.0

##Installation
1. Download the package
2. <code>cd</code> into referenceViz/src/
3. Type <code>make</code>
4. Add the referenceViz/bin/ directory into your PATH

##Run
1. <code>getPlot.py</code>: a wrapper to run various analyses, including: 
    Contiguity 
    Coverage
    N50
    SNP rate
    Indel rate
    Indel length distribution
    CNV
    dbSNP and 1000 Genomes project SNPs and short indels (<= 10bp) validation

    Input include:
    Location of output directory after running [referenceScript pipeline ](https://github.com/benedictpaten/referenceScripts/blob/master/bin/pipeline.py)
    If dbSNP and 1000 Genomes project SNPs validation is included, the dbSNP file with the known SNPs is required. 
    If dbSNP and 1000 Genomes project indels validation is included, the dbSNP file with the known indels is required. 
    The SNP files must be in the tab-separated format of the fields specified [here ](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=varRep&hgta_track=snp135&hgta_table=snp135&hgta_doSchema=describe+table+schema), excluding the "bin" field.

2. <code>mapReadsToRef.py</code>: script to compare mapping of short reads to two different reference sequences (i.e the C. Ref. sequence and GRCh37 sequence)

3. <code>mapLargeIndels.py</code>: script to extract the large-indel sequences and map them to the reference.

    <code>largeIndelTab.py</code>: generates the summary table of the large indel mapping.

4. <code>uniqMapStats.py</code>: script to compare reads that map uniquely to one reference but not uniquely to the other.
   
    <code>uniqMapTab.py</code>: generates the summary table for uniquely-mapping comparisons.

