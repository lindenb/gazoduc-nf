SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20220911.$(BUILD).boxplotdepth
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk

NF=$(realpath boxplotdepth.nf)

all: $(NF) ../../confs/${HOSTNAME}.cfg $(OUTDIR)/bams.txt $(OUTDIR)/jeter.bed $(OUTDIR)/groups.txt
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg  run -lib ../../lib -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--bams "$(OUTDIR)/bams.txt" \
		--bed "$(OUTDIR)/jeter.bed" \
		--sample2group $(OUTDIR)/groups.txt 
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/groups.txt:
	mkdir -p $(dir $@)
	echo -e "S1\tCASE" > $@
	echo -e "S2\tCASE" >> $@
	echo -e "S3\tCASE" >> $@
	echo -e "S4\tCTRL" >> $@
	echo -e "S5\tCTRL" >> $@

$(OUTDIR)/bams.txt : 
	mkdir -p $(dir $@)
	find /LAB-DATA/BiRD/users/lindenbaum-p/src/jvarkit-git/src/test/resources/ -type f -name "S*.bam" > $@	

$(OUTDIR)/jeter.bed:
	mkdir -p $(dir $@)
	echo -e "RF01\t100\t110\tTR1\tEX1" > $@
	echo -e "RF01\t120\t200\tTR1\tEX2" >> $@
	echo -e "RF01\t230\t250\tTR1\tEX4" >> $@
	echo -e "RF03\t100\t1000\tTR2\tEX1" >> $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -lib ../../lib -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
