BUILD=hs37d5
PREFIX=20220624.smoove.
include ../../../data/reference/references.mk
all: work/input.cases work/input.controls ../../confs/login-01.compute.bird2.prive.cfg
	module load nextflow && nextflow -C ../../confs/login-01.compute.bird2.prive.cfg  run -resume delly2.nf \
		--publishDir ${PWD}/work \
		--reference $(REF) \
		--prefix "$(PREFIX)" \
		--cases ${PWD}/work/input.cases \
		--controls ${PWD}/work/input.controls
	-dot -T svg -o workflow.svg work/$(PREFX)workflow.dot

work/input.cases: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

work/input.controls: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep PREGO -m2 $< > $@

README.md: delly2.nf
	module load nextflow && nextflow run $< --help | tail -n+3 > $@

clean:
	rm -rvf work .nextflow .nextflow.l*
