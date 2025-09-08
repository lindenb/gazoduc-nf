## Parameters

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