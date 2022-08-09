SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20220712.$(BUILD).remap
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=$(realpath mosdepth.nf)

all: $(NF) ../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.txt $(OUTDIR)/jeter.bed
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.txt" \
		--bed "$(OUTDIR)/jeter.bed"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.txt : 
	mkdir -p $(dir $@)
	find /LAB-DATA/BiRD/users/lindenbaum-p/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@	

$(OUTDIR)/jeter.bed:
	mkdir -p $(dir $@)
	echo -e "RF01\t100\t1000" > $@
	echo -e "RF03\t100\t1000" >> $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
