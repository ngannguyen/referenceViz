binPath = ../bin
libPath = ../lib

progs = $(notdir $(wildcard indelTable.py indelDistPlot.py getPlots.py contiguityPlot.py coveragePlot.py n50Plot.py snpPlot.py))

.PHONY: %.all

targets = ${progs:%=${binPath}/%} 

all: ${targets}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean :
	rm -f ${progs:%=${binPath}/%}
