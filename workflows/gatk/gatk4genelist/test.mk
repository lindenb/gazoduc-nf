BUILD=hs37d5
PREFIX=20220625.gatk4genes
NF=gatk4genes.nf
OUTDIR=/SCRATCH-BIRD/users/${USER}/work/NEXTFLOW/$(PREFIX)/work
include ../../../../data/reference/references.mk
all: $(OUTDIR)/jeter.bams $(NF) ../../../confs/login-01.compute.bird2.prive.cfg $(OUTDIR)/genes.txt
	module load nextflow && nextflow -C ../../../confs/login-01.compute.bird2.prive.cfg  run -work-dir "$(OUTDIR)" -lib  ../../../lib/  -resume $(NF) \
		--publishDir $(OUTDIR) \
		--reference $(REF) \
		--prefix "$(PREFIX)." \
		--genes "$(OUTDIR)/genes.txt" \
		--mapq 30 \
		--bams $(OUTDIR)/jeter.bams
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/jeter.bams: ../../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

$(OUTDIR)/genes.txt:
	mkdir -p $(dir $@)
	echo "SCN5A,SCN10A" | tr "," "\n" > $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../../confs/login-01.compute.bird2.prive.cfg run $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
