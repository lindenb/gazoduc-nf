BUILD=hs37d5
PREFIX=20220626.gazoduc.test.smoove
OUTDIR=/SCRATCH-BIRD/users/${USER}/work/NEXTFLOW/$(PREFIX)/work
include ../../../data/reference/references.mk

all: $(OUTDIR)/jeter.cases.bams $(OUTDIR)/jeter.ctrls.bams smoove.population.nf ../../confs/${HOSTNAME}.cfg
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" -resume smoove.population.nf \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--cases "$(OUTDIR)/jeter.cases.bams" \
		--controls "$(OUTDIR)/jeter.ctrls.bams" \
		--gff3 $(GFF) \
		--smoove_image /LAB-DATA/BiRD/users/lindenbaum-p/smoove-0.2.6.simg
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/jeter.cases.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

$(OUTDIR)/jeter.ctrls.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E 'PREGO' $< -m2 > $@

README.md: smoove.population.nf
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
