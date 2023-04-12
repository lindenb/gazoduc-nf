SHELL=/bin/bash
BUILD=rotavirus
NF=pggb.01.nf
PREFIX=20230406.$(BUILD).test.pggb
OUTDIR=/SCRATCH-BIRD/users/${USER}/work/NEXTFLOW/$(PREFIX)/work
include ../../../../data/reference/references.mk

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.list $(OUTDIR)/jeter.bed
	module load nextflow/22.04.0 openjdk/11.0.8 && nextflow -C ../../../confs/${HOSTNAME}.cfg  run -lib ../../../lib -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.list" \
		--bed $(OUTDIR)/jeter.bed \
		--conda "${CONDA_ENVS_PATH}"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.list : 
	mkdir -p $(dir $@)
	find $(dir $(REF)) -type f -name "S*.bam" > $@
	
$(OUTDIR)/jeter.bed: $(REF)
	mkdir -p $(dir $@)
	awk -F '\t' '{printf("%s\t0\t%s\n",$$1,$$2);}' "$(REF).fai" > $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
