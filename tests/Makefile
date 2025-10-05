SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GRCh38
.PHONY=all reference test_wgsim test_fastqc test_bwa


all: test_wgsim test_fastqc test_bwa

test_bwa : results/FASTQ/samplesheet.csv base.config reference.config
	which wgsim || true
	echo 'params.wgsim=20' > local.cfg
	$(NEXTFLOW) -c base.config -c reference.config -c local.cfg run ../workflows/bwa/main.nf --samplesheet $<
	rm local.cfg

test_fastqc: results/FASTQ/samplesheet.csv
	$(NEXTFLOW) -c base.config run ../workflows/fastqc/main.nf --samplesheet $<

test_wgsim: results/FASTQ/samplesheet.csv

results/FASTQ/samplesheet.csv: base.config reference.config
	$(NEXTFLOW) -c base.config -c reference.config  run ../workflows/wgsim/main.nf --n_samples 10 --samplesheet $<

reference.config: REF/$(BUILD).fa 
	echo 'params {' > $@
	echo '  fasta="REF/$(BUILD).fa"' >> $@
	echo '  ucsc_name="hg38"' >> $@
	echo '  bwa_index_directory=null' >> $@
	echo '}'  >> $@


reference: REF/$(BUILD).fa REF/$(BUILD).fa.fai REF/$(BUILD).dict

DATA/FASTQ/samplesheet1.csv:
	mkdir -p $(dir $@)
	cd $(dir $@) && \
		wget -O jeter.zip "https://zenodo.org/records/1236641/files/test_fastq_small.zip?download=1"  && \
		unzip jeter.zip && \
		rm jeter.zip
	find $(dir $@) -type f -name "*.fastq" -exec gzip '{}' ';'
	echo "sample,fastq_1,fastq_2" > $@
	find $(dir $@) -name "*.fastq.gz" | sort -V | paste -d, - - | awk '{printf("SAMPLE%d,%s\n",NR,$$0);}' >> $@

REF/$(BUILD).xxfa.bwt:
	mkdir -p $(dir $@)
	cd $(dir $@) && wget -q -O jeter.tar.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_no_alt_analysis_set.fna.bwa_index.tar.gz" && \
		tar xvfz jeter.tar.gz && \
		mv $(BUILD)_no_alt_analysis_set.fna.amb $(BUILD).fa.amb && \
		mv $(BUILD)_no_alt_analysis_set.fna.ann $(BUILD).fa.ann && \
		mv $(BUILD)_no_alt_analysis_set.fna.bwt $(BUILD).fa.bwt && \
		mv $(BUILD)_no_alt_analysis_set.fna.pac $(BUILD).fa.pac && \
		mv $(BUILD)_no_alt_analysis_set.fna.sa $(BUILD).fa.sa


REF/$(BUILD).gff.gz.tbi: REF/$(BUILD).gff.gz
	tabix -p gff -f $<

REF/$(BUILD).gff.gz:
	mkdir -p $(dir $@)
	wget -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz" |\
		gunzip -c | grep '^chr22' | \
		bgzip  > $@


REF/$(BUILD).gtf.gz.tbi: REF/$(BUILD).gtf.gz
	tabix -p gff -f $<

REF/$(BUILD).gtf.gz:
	mkdir -p $(dir $@)
	wget -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz" |\
		gunzip -c | grep '^chr22' | \
		bgzip  > $@



REF/$(BUILD).fa:
	mkdir -p $(dir $@) && \
		wget -O - "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz" |\
		gunzip -c  > $@


