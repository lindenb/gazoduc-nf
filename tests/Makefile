SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GCA_000001405.15_GRCh38
.PHONY=all reference test_fastqc


all: test_fastqc


test_fastqc: DATA/FASTQ/samplesheet1.csv base.config
	$(NEXTFLOW) -c conda.config run ../workflows/fastqc/main.nf --samplesheet $<

reference.config: REF/$(BUILD).fa REF/$(BUILD).fa.fai REF/$(BUILD).dict REF/$(BUILD).fa.bwt
	echo 'params {' > $@
	echo '  fasta="REF/$(BUILD).fa"' >> $@
	echo '  fai="REF/$(BUILD).fa.fai"' >> $@
	echo '  dict="REF/$(BUILD).dict"' >> $@
	echo '  ucsc_name="hg38"' >> $@
	echo '}'  >> $@

base.config:
	echo 'conda.enabled = true'  > $@
	echo 'process {' >> $@
	echo 'maxForks=1' >> $@
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
	find $(dir $@) -name "*.fastq.gz" | sort -V | paste -sd, | awk '{printf("SAMPLE%d,%s\n",NR,$$0);}' >> $@
	
REF/$(BUILD).fa.bwt:
	mkdir -p $(dir $@)
	cd $(dir $@) && wget -O jeter.tar.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_no_alt_analysis_set.fna.bwa_index.tar.gz" && \
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
	wget -O - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_full_analysis_set.refseq_annotation.gff.gz" |\
		gunzip -c |\
		bgzip  > $@


REF/$(BUILD).gtf.gz.tbi: REF/$(BUILD).gtf.gz
	tabix -p gff -f $<
	
REF/$(BUILD).gtf.gz:
	mkdir -p $(dir $@) 
	wget -O - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_full_analysis_set.refseq_annotation.gtf.gz" |\
		gunzip -c |\
		bgzip  > $@

REF/$(BUILD).dict: REF/$(BUILD).fa
        samtools dict -a "$(BUILD)" -A -o $@ $<


REF/$(BUILD).fa.fai: REF/$(BUILD).fa
	mkdir -p $(dir $@) && \
		wget -O - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_no_alt_analysis_set.fna.fai" > $@


REF/$(BUILD).fa:
	mkdir -p $(dir $@) && \
		wget -O - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$(BUILD)/seqs_for_alignment_pipelines.ucsc_ids/$(BUILD)_no_alt_analysis_set.fna.gz" |\
		gunzip -c  > $@

install.conda:


	touch $@


