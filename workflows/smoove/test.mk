BUILD=hs37d5
PREFIX=202206.smoove.
include ../../../data/reference/references.mk
all: work/jeter.cases.bams work/jeter.ctrls.bams smoove.population.nf ../../confs/login-01.compute.bird2.prive.cfg
	module load nextflow && nextflow -C ../../confs/login-01.compute.bird2.prive.cfg  run -resume smoove.population.nf \
		--publishDir ${PWD}/work \
		--reference $(REF) \
		--prefix $(PREFIX) \
		--cases ${PWD}/work/jeter.cases.bams \
		--controls ${PWD}/work/jeter.ctrls.bams \
		--gff3 $(GFF) \
		--smoove_image /LAB-DATA/BiRD/users/lindenbaum-p/smoove-0.2.6.simg
	-dot -T svg -o workflow.svg work/$(PREFIX)workflow.dot

work/jeter.cases.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E '(B00ICUM|B00HOR6)' $< > $@

work/jeter.ctrls.bams: ../../../2021/20210526.resources/Resources.hs37d5.wgs.bams.list
	mkdir -p $(dir $@)
	grep -E 'PREGO' $< -m2 > $@

README.md: smoove.population.nf
	module load nextflow && nextflow -C ../../confs/login-01.compute.bird2.prive.cfg  run $< --help | tail -n+3 > $@

clean:
	rm -rvf work .nextflow .nextflow.*
