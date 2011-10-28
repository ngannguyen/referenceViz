binPath = ../bin
libPath = ../lib

progs = $(notdir $(wildcard indelCheckPlot.py snpCheckPlot.py indelPlot.py overlapSnps.py mergeDbsnpsAnd1ksnps.py mappingPlot.py bamGetPairs.py indelStats.py snpStats.py parseBam.py getReads.py mapReadsToRef.py cnvPlot.py summaryStats.py indelTable.py indelDistPlot.py getPlots.py contiguityPlot.py coveragePlot.py n50Plot.py snpPlot.py))

.PHONY: %.all

targets = ${progs:%=${binPath}/%} 

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
