SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20221121.$(BUILD).plot
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=wgscovplot.01.nf


all: $(NF) ../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.txt
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.txt" \
		--maxcov 10 \
		--include_contig_regex ".*"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.txt : 
	mkdir -p $(dir $@)
	find /LAB-DATA/BiRD/users/lindenbaum-p/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@	


README.md: $(NF)
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
