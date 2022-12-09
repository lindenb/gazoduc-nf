SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20220712.$(BUILD).remap
OUTDIR=$(PWD)/work
include ../../../../data/reference/references.mk

NF=$(realpath graphtyper.genotype.nf)

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.list $(OUTDIR)/beds.list
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--depth_method samtoolsdepth \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.list" \
		--beds "$(OUTDIR)/beds.list"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.list : 
	mkdir -p $(dir $@)
	find ${HOME}/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@	

$(OUTDIR)/beds.list: $(OUTDIR)/jeter.bed
	echo $< > $@

$(OUTDIR)/jeter.bed : $(addsuffix .fai, $(REF))
	mkdir -p $(dir $@)
	awk '{printf("%s\t0\t%d\n",$$1,$$2);}' $< > $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
