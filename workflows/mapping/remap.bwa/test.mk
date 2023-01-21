SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20220712.$(BUILD).remap
OUTDIR=$(PWD)/work
include ../../../../data/reference/references.mk

NF=$(realpath remap.bwa.nf)

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.list
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference_in $(REF) \
		--reference_out $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.list"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.list :
	mkdir -p $(dir $@)
	find ${HOME}/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@


README.md: $(NF)
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg -C remap.bwa.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
