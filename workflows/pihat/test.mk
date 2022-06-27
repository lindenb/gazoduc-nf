SHELL=/bin/bash
BUILD=hs37d5
VCF=/LAB-DATA/BiRD/shares/ITX/u1087/lindenb/20210624.hs37d5.NTS271/20210625.hs37d5.NTS271.merged.vcf.gz
PREFIX=20220626.gazoduc.test.pihat
OUTDIR=/SCRATCH-BIRD/users/${USER}/work/NEXTFLOW/$(PREFIX)/work
include ../../../data/reference/references.mk

all: pihat.nf ../../confs/${HOSTNAME}.cfg $(OUTDIR)/samples.txt $(VCF)
	$(MAKE) -C ../../src
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg  run -lib ../../lib -work-dir "$(OUTDIR)" -resume pihat.nf \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--vcf "$(VCF)" \
		--samples $(OUTDIR)/samples.txt
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/samples.txt: $(VCF)
	mkdir -p $(dir $@)
	bcftools query -l $< > $@

README.md: pihat.nf
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
