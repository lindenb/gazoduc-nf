
## About

call MT genome

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) path to multiple bams files [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow \
	--publishDir output \
	--prefix "analysis." \
	--reference /path/to/reference.fasta \
	--vcf in.vcf.gz \
	--bed /path/to/in.bed
```

## Workflow

![workflow](./workflow.svg)
  
## See also


