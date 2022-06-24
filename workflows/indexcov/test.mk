BUILD=hs37d5
PREFIX=202206.indexcov.
include ../../../data/reference/references.mk
all: work/jeter.bams indexcov.nf ../../confs/login-01.compute.bird2.prive.cfg
	module load nextflow && nextflow -C ../../confs/login-01.compute.bird2.prive.cfg  run -resume indexcov.nf \
		--publishDir ${PWD}/work \
		--reference $(REF) \
		--prefix $(PREFIX) \
		--mapq 30 \
		--bams ${PWD}/work/jeter.bams
	-dot -T svg -o workflow.svg work/$(PREFIX)workflow.dot

work/jeter.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

README.md: indexcov.nf
	module load nextflow && nextflow run $< --help | tail -n+3 > $@

clean:
	rm -rvf work .nextflow .nextflow.*
