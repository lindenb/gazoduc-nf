SHELL=/bin/bash
BUILD=hs37d5
PREFIX=20220731.$(BUILD).retrocopies
OUTDIR=$(PWD)/work
include ../../../../data/reference/references.mk

VCF=${HOME}/externalshares/GENMED/20220507.genmed.hs37d5.smoove.all.bcf

NF=$(realpath vcf.retrocopy.nf)

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(VCF) $(GTF) $(REF)
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--gtf "$(GTF)"  \
		--vcf $(VCF)
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

README.md: $(NF)
	module load nextflow && nextflow -C ../../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
