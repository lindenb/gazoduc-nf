SHELL=/bin/bash
BUILD=hs38me
PREFIX=20221205.$(BUILD).cardiobed
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=cardiobed.nf

all: $(NF) ../../confs/${HOSTNAME}.cfg
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--gtf $(GFF)
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
