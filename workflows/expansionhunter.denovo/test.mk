BUILD=hs37d5
PREFIX=20221017.xhunter
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=expansionhunter.de.novo.01.nf
all: $(NF) $(OUTDIR)/jeter.bams ../../confs/${HOSTNAME}.cfg
	module load java/1.8.0_131 nextflow/22.04.0  && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--cases "$(OUTDIR)/cases.bams"
		--controls "$(OUTDIR)/controls.bams"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/cases.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E 'PREGO' $< | head -n 10  > $@


README.md: $(NF)
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
