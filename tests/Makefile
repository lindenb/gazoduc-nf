SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GRCh38
.PHONY=all hard_tests reference test_wgsim test_fastqc test_bwa test_hapcaller stub

all: stub

stub: data/rotavirus_rf.fa data/rotavirus_rf.vcf.gz
	mkdir -p STUB
	$(NEXTFLOW) run ../workflows/wgsim/main.nf  -work-dir STUB/WGSIM  -resume -stub --outdir STUB/WGSIM/results --fasta $<
	$(NEXTFLOW) run ../workflows/fastqc/main.nf  -work-dir STUB/FASTQC  -resume -stub --outdir STUB/FASTQC/results  --samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv
	$(NEXTFLOW) run ../workflows/bwa/main.nf  -work-dir STUB/BWA  -resume -stub --outdir STUB/BWA/results  \
		--samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv \
		--fasta $< \
		--known_sites data/rotavirus_rf.vcf.gz
	# combinegvcfs
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $< \
		--dbsnp data/rotavirus_rf.vcf.gz
	# glnexus
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $< \
		--gvcf_merge_method glnexus \
		--dbsnp data/rotavirus_rf.vcf.gz
	# method bam2vcf
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $< \
		--method bam2vcf \
		--dbsnp data/rotavirus_rf.vcf.gz
	# method direct
	$(NEXTFLOW) run ../workflows/gatk/gatk4gvcfs3/main.nf  -work-dir STUB/HC  -resume -stub --outdir STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $< \
		--method direct \
		--dbsnp data/rotavirus_rf.vcf.gz
	# somalier
	$(NEXTFLOW) run ../workflows/somalier/somalier.bams/main.nf  -work-dir STUB/SOMALIER -resume -stub --outdir STUB/SOMALIER/results \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta data/rotavirus_rf.fa
	# sashimi encode
	touch jeter.gtf.gz jeter.gtf.gz.tbi
	echo -e 'RF01\t1\t1000' > jeter.bed
	$(NEXTFLOW) run ../workflows/sashimi.encode/main.nf  -work-dir STUB/ENCODE -stub --outdir STUB/ENCODE \
		--bed jeter.bed \
		--gtf jeter.gtf.gz \
		--fasta data/rotavirus_rf.fa \
		--ucsc_name "hg38"
	# delly
	$(NEXTFLOW) run ../workflows/delly2/main.nf  -work-dir STUB/DELLY -stub --outdir STUB/DELLY \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta data/rotavirus_rf.fa
	# downsample
	$(NEXTFLOW) run ../workflows/downsample.simulation/main.nf -work-dir STUB/LOWPASS -stub -resume --outdir STUB/LOWPASS \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv  \
		--fasta data/rotavirus_rf.fa

	
hard_tests: test_wgsim test_fastqc test_bwa test_hapcaller stub

test_hapcaller: work/BWA/results/BAMS/samplesheet.csv
	$(NEXTFLOW) run -work-dir work/HC -profile gitlabci -c reference.config  -c nproc.cfg  ../workflows/gatk/gatk4gvcfs3/main.nf --outdir work/HC/results --samplesheet $<
	find work/HC/ -type f


test_bwa : work/BWA/results/BAMS/samplesheet.csv


work/BWA/results/BAMS/samplesheet.csv : work/WGSIM/results/FASTQ/samplesheet.csv  reference.config nproc.cfg REF/$(BUILD).dbsnp.vcf.gz
	$(NEXTFLOW) run -work-dir work/BWA -profile gitlabci -c reference.config  -c nproc.cfg  ../workflows/bwa/main.nf --outdir work/BWA/results --samplesheet $< --known_sites REF/$(BUILD).dbsnp.vcf.gz
	find work/BWA/results -type f

test_fastqc: work/WGSIM/results/FASTQ/samplesheet.csv nproc.cfg
	$(NEXTFLOW) run -work-dir work/FASTQC -profile gitlabci -c nproc.cfg  --outdir FASTQC ../workflows/fastqc/main.nf --samplesheet $<


test_wgsim: work/WGSIM/results/FASTQ/samplesheet.csv

work/WGSIM/results/FASTQ/samplesheet.csv: reference.config nproc.cfg
	$(NEXTFLOW) run  -work-dir work/WGSIM -profile gitlabci  -c reference.config  -c nproc.cfg   --outdir work/WGSIM/results  ../workflows/wgsim/main.nf --n_samples 15 --samplesheet $<
	
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
