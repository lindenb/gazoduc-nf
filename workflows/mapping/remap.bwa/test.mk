SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20230121.$(BUILD).remap
OUTDIR=$(PWD)/work
include ../../../../data/reference/references.mk

NF=$(realpath remap.bwa.nf)

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.list
	module load nextflow/22.04.0 openjdk/11.0.8 && nextflow -C ../../../confs/${HOSTNAME}.cfg run -lib ../../../lib -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference_in $(REF) \
		--reference_out $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.list" \
		--split_fastqs_count 2 \
		--split_fastq_ignore_if_size 10
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/bams.list :
	mkdir -p $(dir $@)
	find /LAB-DATA/BiRD/users/lindenbaum-p/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@


README.md: $(NF)
	module load nextflow/22.04.0 && nextflow -C ../../../confs/${HOSTNAME}.cfg  run -lib ../../../lib -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.* trace.tx*
