SHELL=/bin/bash
BUILD=rotavirus
PREFIX=20220712.$(BUILD).remap
OUTDIR=$(PWD)/work
include ../../../../data/reference/references.mk

NF=$(realpath fastmap.01.nf)

all: $(NF) ../../../confs/${HOSTNAME}.cfg $(OUTDIR)/input.table
	module load nextflow/22.04.0 && nextflow -C ../../../confs/${HOSTNAME}.cfg  run -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--reference $(REF) \
		--prefix "$(PREFIX)."  \
		--fastqs "$(OUTDIR)/input.table"
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

$(OUTDIR)/input.table: 
	mkdir -p $(dir $@)
	echo -e "sample\tR1\tR2" > $@
	find ${HOME}/src/jvarkit-git/src/test/resources/ -type f -name "S*fq.gz" |  sort | paste - - | awk '{printf("S%s\t%s\n",NR,$$0);}' >> $@
	

README.md: $(NF)
	module load nextflow/22.04.0 && nextflow -C ../../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
