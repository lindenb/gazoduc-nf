BUILD=rotavirus
SHELL=/bin/bash
PREFIX=20230103.$(BUILD).etching
OUTDIR=$(PWD)/work
include ../../../data/reference/references.mk
NF=etching.germline.nf

RSRC=$(dir $(REF))

all: $(NF) ../../confs/${HOSTNAME}.cfg $(OUTDIR)/fastqs.csv
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--prefix "$(PREFIX)."  \
		--reference "$(REF)"  \
		--fastqs "$(OUTDIR)/fastqs.csv"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/fastqs.csv:
	mkdir -p $(dir $@)
	echo "sample,R1,R2" > $@
	echo "S1,$(RSRC)/S1.R1.fq.gz,$(RSRC)/S1.R2.fq.gz" >> $@
	echo "S2,$(RSRC)/S2.R1.fq.gz,$(RSRC)/S2.R2.fq.gz" >> $@
	echo "S3,$(RSRC)/S3.R1.fq.gz,$(RSRC)/S3.R2.fq.gz" >> $@

README.md: $(NF)
	module load nextflow && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
