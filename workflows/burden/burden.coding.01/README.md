
## About

Burden for coding regions.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --cases <file> file containing the cases' names. One per line.
  * --controls <file> file containing the controls' names. One per line.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume burden.coding.01.nf \
        --publishDir output \
        --prefix "analysis." \
        --reference /path/to/reference.fasta \
        --vcf /path/to/my.vcf.gz \
        --cases /path/to/cases.file \
        --controls /path/to/controls.file
```

## Workflow

![workflow](./workflow.svg)


