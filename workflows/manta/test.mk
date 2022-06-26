BUILD=hs37d5
PREFIX=20220626.gazoduc.test.manta
OUTDIR=/SCRATCH-BIRD/users/${USER}/work/NEXTFLOW/$(PREFIX)/work
include ../../../data/reference/references.mk
all: $(OUTDIR)/bams.list
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" -resume manta.nf \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--bams "$(OUTDIR)/bams.list"
	 -dot -T svg -o workflow.svg "$(OUTDIR)/$(PREFIX).workflow.dot"

$(OUTDIR)/bams.list : ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

README.md: manta.nf
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.lo*
