SHELL=/bin/bash
NEXTFLOW?=nextflow
BUILD=GRCh38
PROFILE?= -profile gitlabci
.PHONY=all hard_tests reference test_wgsim test_fastqc test_bwa test_hapcaller test_downsample stub nf-test

all: work/BWA/results/BAMS/samplesheet.csv
	nf-test test --profile micromamba --verbose

nf-test: nf-test-wgsim

nf-test-wgsim:
	cd .. && nf-test test --profile micromamba tests/workflows/wgsim/main.nf.test


stub: data/rotavirus_rf.fa data/rotavirus_rf.vcf.gz
	mkdir -p STUB
	$(NEXTFLOW) run ../workflows/wgsim/main.nf  -work-dir STUB/WGSIM  -resume -stub --outdir STUB/WGSIM/results --fasta $<
	$(NEXTFLOW) run ../workflows/fastqc/main.nf  -work-dir STUB/FASTQC  -resume -stub --outdir STUB/FASTQC/results  --samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv
	$(NEXTFLOW) run ../workflows/bwa/main.nf  -work-dir STUB/BWA  -resume -stub --outdir STUB/BWA/results  \
		--samplesheet STUB/WGSIM/results/FASTQ/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--known_sites  $${PWD}/data/rotavirus_rf.vcf.gz
	# pihat
	$(NEXTFLOW) run ../workflows/pihat/main.nf  -work-dir STUB/PIHAT  -resume -stub --outdir STUB/PIHAT/results  \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--vcf  $${PWD}/data/rotavirus_rf.vcf.gz	
	# combinegvcfs
	$(NEXTFLOW) run ../workflows/gatk/haplotypecaller/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--dbsnp $${PWD}/data/rotavirus_rf.vcf.gz
	# glnexus
	$(NEXTFLOW) run ../workflows/gatk/haplotypecaller/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--gvcf_merge_method glnexus \
		--dbsnp data/rotavirus_rf.vcf.gz
	# method bam2vcf
	$(NEXTFLOW) run ../workflows/gatk/haplotypecaller/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--method bam2vcf \
		--dbsnp  $${PWD}/data/rotavirus_rf.vcf.gz
	# method bam2vcf dnc
	$(NEXTFLOW) run ../workflows/gatk/haplotypecaller/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
		--samplesheet STUB/BWA/results/BAMS/samplesheet.csv \
		--fasta $${PWD}/data/rotavirus_rf.fa \
		--method bam2vcf \
		--divide_and_conquer true \
		--dbsnp  $${PWD}/data/rotavirus_rf.vcf.gz
	# method direct
	$(NEXTFLOW) run ../workflows/gatk/haplotypecaller/main.nf  -work-dir STUB/HC  -resume -stub --outdir $${PWD}/STUB/HC/results  \
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
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/haplotypecaller/main.nf \
		--outdir $${PWD}/work/HC/results \
		--samplesheet $<
	find work/HC/results -type f
	#
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/haplotypecaller/main.nf \
		--outdir $${PWD}/work/HC/results \
		--gvcf_merge_method glnexus \
		--samplesheet $<
	find work/HC/results -type f
	#
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/haplotypecaller/main.nf \
		--outdir $${PWD}/work/HC/results \
		--method direct \
		--samplesheet $<
	find work/HC/results -type f
	#
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/haplotypecaller/main.nf \
		--outdir $${PWD}/work/HC/results \
		--method bam2vcf \
		--samplesheet $<
	find work/HC/results -type f
	#
	$(NEXTFLOW) run -resume -work-dir work/HC $(PROFILE) -c reference.config  -c nproc.cfg  ../workflows/gatk/haplotypecaller/main.nf \
		--outdir $${PWD}/work/HC/results \
		--method bam2vcf --gvcf_merge_method glnexus \
		--samplesheet $<
	find work/HC/results -type f



test_bwa : work/BWA/results/BAMS/samplesheet.csv



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


reference: REF/$(BUILD).fa REF/$(BUILD).fa.fai REF/$(BUILD).dict

	
REF/$(BUILD)/BAMS/samplesheet.csv : REF/$(BUILD)/FASTQ/samplesheet.csv REF/$(BUILD).dbsnp.vcf.gz
	$(NEXTFLOW) run -resume \
		-work-dir work/BWA \
		-c REF/$(BUILD).config \
		-profile gitlabci,$(BUILD)CI \
		../workflows/bwa/main.nf \
		--outdir $${PWD}/REF/$(BUILD) \
		--samplesheet $< \
		--known_sites $${PWD}/REF/$(BUILD).dbsnp.vcf.gz
	find  $${PWD}/REF/$(BUILD) -type f -name "*.am" || true
	find  $${PWD}/REF/$(BUILD) -type f -name "*.csv" || true
	cut -d, -f 1-8 $@ > $(addsuffix .tmp,$@) && mv $(addsuffix .tmp,$@) $@
	cat $@


REF/$(BUILD)/FASTQ/samplesheet.csv: REF/$(BUILD).config
	$(NEXTFLOW) run -resume \
		-work-dir work/WGSIM \
		../workflows/wgsim/main.nf \
		-c $<  \
		-profile gitlabci,$(BUILD)CI \
		--outdir $${PWD}/REF/$(BUILD) \
		--n_samples 7
	find  $(dir $@) -type f
	cat $@


REF/$(BUILD).config : REF/$(BUILD)fa.bwt REF/$(BUILD).gff.gz.tbi REF/$(BUILD).gtf.gz.tbi REF/$(BUILD).fa REF/nproc.cfg
	mkdir -p $(dir $@)
	echo "profiles {" > $@
	echo "${BUILD}CI {" >> $@
	echo "params {" >> $@
	echo "  fasta='$${PWD}/REF/$(BUILD).fa'" >> $@
	echo '  ucsc_name="hg38"' >> $@
	echo "  bwa_index_directory='$${PWD}/REF'" >> $@
	echo '  gff="$(PWD)/REF/$(BUILD).gff.gz"' >> $@
	echo "}" >> $@
	echo "}" >> $@
	echo "}" >> $@
	echo 'includeConfig("nproc.cfg")' >> $@
	cat $@

REF/nproc.cfg:
	mkdir -p $(dir $@)
	nproc --all | awk '{N=int($$0);printf("process {\nmaxForks=%d\n}\n",(N<3?1:N-3));}' > $@
	cat $@


REF/$(BUILD)fa.bwt: REF/$(BUILD).fa
	mkdir -p $(dir $@)
	bwa index $<

REF/$(BUILD).gff.gz.tbi: REF/$(BUILD).gff.gz
	tabix -p gff -f $<

REF/$(BUILD).gff.gz:
	mkdir -p $(dir $@)
	wget -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz" |\
		gunzip -c | grep '^chr22' | \
		LC_ALL=C sort -t $$'\t' -k1,1 -k4,4n -T . |\
		bgzip  > $@


REF/$(BUILD).gtf.gz.tbi: REF/$(BUILD).gtf.gz
	tabix -p gff -f $<

REF/$(BUILD).gtf.gz:
	mkdir -p $(dir $@)
	wget -qO - "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz" |\
		gunzip -c | grep '^chr22' | \
		LC_ALL=C sort -t $$'\t' -k1,1 -k4,4n -T . |\
		bgzip  > $@

REF/$(BUILD).fa:
	mkdir -p $(dir $@) && \
		wget --progress=dot:mega -O $(addsuffix .gz,$@) "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz"
	gunzip -f  $(addsuffix .gz,$@)
	touch -c $@

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
