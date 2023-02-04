BUILD=rotavirus
PREFIX=20220626.gazoduc.test.manta
OUTDIR=${PWD}/work
include ../../../data/reference/references.mk
all: $(OUTDIR)/bams.list
	module load nextflow/22.04.0 openjdk/11.0.8 && nextflow -C ../../confs/${HOSTNAME}.cfg run -lib ../../lib -work-dir "$(OUTDIR)" -resume manta.nf \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--bams "$(OUTDIR)/bams.list"
	 -dot -T svg -o workflow.svg "$(OUTDIR)/$(PREFIX).workflow.dot"

$(OUTDIR)/bams.list : 
	mkdir -p $(dir $@)
	 find $(realpath ../../../../src/jvarkit/src/test/resources/) -type f -name "S*.bam" > $@

README.md: manta.nf
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.lo*
