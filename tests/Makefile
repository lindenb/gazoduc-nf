SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GRCh38
.PHONY=all reference test_wgsim test_fastqc test_bwa


all: test_wgsim test_fastqc test_bwa

test_bwa : DATA/FASTQ/samplesheet2.csv

DATA/FASTQ/samplesheet2.csv : results/FASTQ/samplesheet.csv  reference.config nproc.cfg REF/$(BUILD).dbsnp.vcf.gz
	$(NEXTFLOW) -profile gitlabci -c reference.config  -c nproc.cfg run ../workflows/bwa/main.nf --outdir BWA --samplesheet $< --known_sites REF/$(BUILD).dbsnp.vcf.gz
	
test_fastqc: results/FASTQ/samplesheet.csv nproc.cfg

	$(NEXTFLOW) -profile gitlabci -c nproc.cfg run ../workflows/fastqc/main.nf --samplesheet $<

	
test_wgsim: results/FASTQ/samplesheet.csv

results/FASTQ/samplesheet.csv: reference.config nproc.cfg
	echo 'params.wgsim=20' > local.cfg
	$(NEXTFLOW) -profile gitlabci  -c reference.config  -c nproc.cfg run ../workflows/wgsim/main.nf --n_samples 10 --samplesheet $<
	rm local.cfg
	
reference.config: REF/$(BUILD).fa 
	echo 'params {' > $@
	echo '  fasta="REF/$(BUILD).fa"' >> $@
	echo '  ucsc_name="hg38"' >> $@
	echo '  bwa_index_directory=null' >> $@
	echo '}'  >> $@

nproc.cfg:
	nproc --all | awk '{N=int($$0);printf("process {\nmaxForks=%d\n}\n",(N<3?1:N-3));}' > $@
	cat $@


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
		wget -qO - "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz" |\
		gunzip -c  > $@

REF/$(BUILD).dbsnp.vcf.gz:
	mkdir -p $(dir $@) && \
	wget -qO - "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" |\
		gunzip -c |\
		grep -v '^##contig=' |\
		head -n 10000 |\
		awk '/^#CHROM/ {printf("##contig=<ID=chr22,length=50818468>\n");}{print;}'|\
		sed 's/^22/chr22/' | bcftools view -O z -o $@ && bcftools index -f -t $@

