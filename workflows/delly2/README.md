
## About

Detects CNV/SV using delly2

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France.

## Options

  * --reference (fasta) indexed fasta reference [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --mapq (int)  min mapping quality . If it's lower than 0 (this is the default) just use the bam index as is. Otherwise, rebuild the bai
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --goleft_version (string) default: "v0.2.4"

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume indexcov.nf \
	--publishDir output \
	--prefix "analysis." \
	--reference /path/to/reference.fasta \
	--bams /path/to/bams.list \
	--mapq 30
```
  
## See also



