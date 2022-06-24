
## About

Detects CNV/SV using manta.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM. [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --survivor_merge_params (string) or empty to disable survivor. [1000 2 1 1 0 30]

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume manta.nf \
	--publishDir output \
	--prefix "analysis." \
	--reference /path/to/reference.fasta \
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

  * https://github.com/Illumina/manta


