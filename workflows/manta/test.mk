BUILD=hs37d5
PREFIX=20220624.manta.
include ../../../data/reference/references.mk
all: work/bams.list
	module load nextflow && nextflow -C ../../confs/login-01.compute.bird2.prive.cfg  run -resume manta.nf \
		--publishDir ${PWD}/work \
		--reference $(REF) \
		--prefix "$(PREFIX)" \
		--bams ${PWD}/work/bams.list
	 -dot -T svg -o workflow.svg "work/$(PREFIX)workflow.dot"

work/bams.list : ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

README.md: manta.nf
	module load nextflow && nextflow run $< --help | tail -n+3 > $@

clean:
	rm -rvf work .nextflow .nextflow.lo*
