SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GRCh38
PROFILE?= -profile gitlabci
.PHONY=all hard_tests reference test_wgsim test_fastqc test_bwa test_hapcaller test_downsample stub

all: stub

stub: data/rotavirus_rf.fa data/rotavirus_rf.vcf.gz
	mkdir -p STUB
	$(NEXTFLOW) run ../workflows/wgsim/main.nf  -work-dir STUB/WGSIM  -resume -stub --outdir STUB/WGSIM/results --fasta $<
	$(NEXTFLOW) run ../workflows/fastqc/main.nf  -work-dir STUB/FASTQC  -resume -stub --outdir STUB/FASTQC/results  --samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv
	$(NEXTFLOW) run ../workflows/bwa/main.nf  -work-dir STUB/BWA  -resume -stub --outdir STUB/BWA/results  \
		--samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--known_sites  $${PWD}/data/rotavirus_rf.vcf.gz
	# combinegvcfs
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--dbsnp $${PWD}/data/rotavirus_rf.vcf.gz
	# glnexus
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--gvcf_merge_method glnexus \
		--dbsnp data/rotavirus_rf.vcf.gz
	# method bam2vcf
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--method bam2vcf \
		--dbsnp  $${PWD}/data/rotavirus_rf.vcf.gz
	# method direct
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--method direct \
		--dbsnp $${PWD}/data/rotavirus_rf.vcf.gz
	# somalier
	$(NEXTFLOW) run ../workflows/somalier/somalier.bams/main.nf  -work-dir STUB/SOMALIER -resume -stub --outdir $${PWD}/STUB/SOMALIER/results \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta $${PWD}/data/rotavirus_rf.fa
	# sashimi encode
	touch jeter.gtf.gz jeter.gtf.gz.tbi
	echo -e 'RF01\t1\t1000' > jeter.bed
	$(NEXTFLOW) run ../workflows/sashimi.encode/main.nf  -work-dir STUB/ENCODE -stub --outdir $${PWD}/STUB/ENCODE \
		--bed jeter.bed \
		--gtf jeter.gtf.gz \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--ucsc_name "hg38"
	# delly
	$(NEXTFLOW) run ../workflows/delly2/main.nf  -work-dir STUB/DELLY -stub --outdir STUB/DELLY \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta $${PWD}/data/rotavirus_rf.fa
	# downsample
	$(NEXTFLOW) run ../workflows/downsample.simulation/main.nf -work-dir STUB/LOWPASS -stub -resume --outdir $${PWD}/STUB/LOWPASS \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta $${PWD}/data/rotavirus_rf.fa

	
hard_tests: test_wgsim test_fastqc test_bwa test_hapcaller test_downsample stub 


test_downsample: work/BWA/results/BAMS/samplesheet.csv
	$(NEXTFLOW) run -resume -work-dir work/LOWPASS $(PROFILE) -c reference.config  -c nproc.cfg  \
		../workflows/downsample.simulation/main.nf \
		--outdir $${PWD}/work/LOWPASS/results --samplesheet $<

test_hapcaller: work/BWA/results/BAMS/samplesheet.csv
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/gatk4gvcfs3/main.nf --outdir $${PWD}/work/HC/results --samplesheet $<
	find work/HC/ -type f


test_bwa : work/BWA/results/BAMS/samplesheet.csv


work/BWA/results/BAMS/samplesheet.csv : work/WGSIM/results/FASTQ/samplesheet.csv  reference.config nproc.cfg REF/$(BUILD).dbsnp.vcf.gz
	$(NEXTFLOW) run -resume -work-dir work/BWA $(PROFILE) -c reference.config  -c nproc.cfg \
		../workflows/bwa/main.nf \
		--outdir $${PWD}/work/BWA/results \
		--samplesheet $< --known_sites $${PWD}/REF/$(BUILD).dbsnp.vcf.gz
	find work/BWA/results -type f
	cut -d, -f 1-8 $@ > $(addsuffix .tmp,$@) && mv $(addsuffix .tmp,$@) $@
	cat $@

test_fastqc: work/WGSIM/results/FASTQ/samplesheet.csv nproc.cfg
	$(NEXTFLOW) run -resume -work-dir work/FASTQC $(PROFILE) -c nproc.cfg  \
		--outdir $${PWD}/work/FASTQC ../workflows/fastqc/main.nf \
		--samplesheet $<


test_wgsim: work/WGSIM/results/FASTQ/samplesheet.csv

work/WGSIM/results/FASTQ/samplesheet.csv: reference.config nproc.cfg
	$(NEXTFLOW) run -resume -work-dir work/WGSIM $(PROFILE) -c reference.config  -c nproc.cfg   \
		--outdir $${PWD}/work/WGSIM/results  ../workflows/wgsim/main.nf --n_samples 7
	cat $@
	
reference.config: REF/$(BUILD).fa 
	echo 'params {' > $@
	echo "  fasta='$${PWD}/REF/$(BUILD).fa'" >> $@
	echo '  ucsc_name="hg38"' >> $@
	echo '  bwa_index_directory=null' >> $@
	echo '}'  >> $@
	cat $@

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


REF/rotavirus_rf.fa:
	mkdir -p $(dir $@)
	

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

clean:
	rm -rf STUB
