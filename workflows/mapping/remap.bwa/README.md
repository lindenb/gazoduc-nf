
## About

Remap bam/cram using bwa.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France

## Options

  * --reference_in (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --reference_out (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume remap.bwa.nf \
	--publishDir output \
	--prefix "analysis." \
	--reference_in /path/to/referenceX.fasta \
	--reference_in /path/to/referenceY.fasta \
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also



