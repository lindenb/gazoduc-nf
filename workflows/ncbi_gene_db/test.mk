SHELL=/bin/bash
PREFIX=20230509.ncbi.genes
OUTDIR=$(PWD)/work

NF=build_ncbi_gene_db.nf


all: $(NF) ../../confs/${HOSTNAME}.cfg
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg  run  -lib ../../lib -work-dir "$(OUTDIR)" -resume $(NF) \
		--publishDir "$(OUTDIR)" \
		--prefix "$(PREFIX)."
	-dot -T svg -o workflow.svg $(OUTDIR)/$(PREFIX).workflow.dot

README.md: $(NF)
	module load nextflow/22.04.0 && nextflow -C ../../confs/${HOSTNAME}.cfg run -work-dir "$(OUTDIR)" $< --help | tail -n+3 > $@

clean:
	rm -rvf "$(OUTDIR)" .nextflow .nextflow.*
