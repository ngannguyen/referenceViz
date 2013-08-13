binPath = ../bin
libPath = ../lib

progs = $(notdir $(wildcard halGeneCompare.py geneCompare.py bedComparator.py makeAssemblyHubFiles.py snpPlot_all2all.py crefKmer.py pslUnion.py geneStats_aggregate.py geneStats.py geneListsIntersect.py gff2psl.py faRename.py vcf1k_hom2het.py falsePosSnps.py mapGtfExons.py vcfFilter.py vcf2constraints.py cigarConvertRef.py gtf2align.py queryStats.py mapLargeIndelsStats.py falsePosSNPsTable.py largeIndelTab.py getSeqRepeatFraction.py uniqMapTab.py bam_convertCoor.py uniqMapStats.py mapLargeIndels.py nonLinearBp.py bam_cmpMapping.py getBamUniqMappedReads.py indelCheckPlot.py snpCheckPlot.py indelPlot.py overlapSnps.py mergeDbsnpsAnd1ksnps.py mappingPlot.py bamGetPairs.py indelStats.py snpStats.py parseBam.py getReads.py mapReadsToRef.py cnvPlot.py summaryStats.py indelTable.py indelDistPlot.py getPlots.py contiguityPlot.py coveragePlot.py n50Plot.py snpPlot.py))

.PHONY: %.all

targets = ${progs:%=${binPath}/%} 

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
