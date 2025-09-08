## Parameters
* `--fasta` : required path to fasta file
* `--fai` : required path to fasta.fai file
* `--dict` : required path to fasta.dict file
* `--ucsc_name` : UCSC name of the genome (hg19/hg38/...)
* `--gtf` : optional path to a tabix indexed GTF file. Will be used for statistics
* `--outdir` : output directory. Default is "results"
* `--hts_type`: type of analysis. Possible values: "WES|WGS"; Default is "WGS". If WES, some SV analysis like DELLY will be disabled.

* `--bed` . Optional. Limit analysis to that BED file. If not defined, use the Whole genome as BED. Chromosome name must match the chromosomes in the reference file.


* `--sampleinfo` ` . Optional. Extra meta data about the samples that could be added to the sample sheet (sex, status, etc...). Those information might be used for example in the QC to show the difference between male/female, case/controls, etc... Samples must be unique. It's CSV file with the following required column in the header : 'sample'

```
sample,sex,status
CD12345,male,case
CD12346,female,control
```

* `--extra_bams` . Optional. And those extra BAM/CRAM files that will be added to the calling/QC. It's CSV file with the following required columns in the header : sample,bam,bai . Files MUST be mapped on the very same FASTA reference. Samples must be unique.

```
sample,bam,bai
CD12345,/path/to/CD12345.bam,/path/to/CD12345.bam.bai
CD12346,/path/to/CD12346.cram,/path/to/CD12346.cram.crai
```