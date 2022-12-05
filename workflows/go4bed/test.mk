SHELL=/bin/bash
BUILD=hs38me
PREFIX=20221205.$(BUILD).go4bed
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=go4bed.nf

all: $(NF) ../../confs/${HOSTNAME}.cfg
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--goTerms "GO:0034765,GO:0043269,GO:0005216,GO:0006811,GO:0034220,GO:0008016" \
		--excludeGoTerms "GO:0045202,GO:0099536,GO:0060078,GO:0003014" \
		--gtf $(GFF)
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
