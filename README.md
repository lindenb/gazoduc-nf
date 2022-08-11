GAZODUC
=======

![Last commit](https://img.shields.io/github/last-commit/lindenb/gazoduc-nf.png)
![works on my machine](https://img.shields.io/static/v1?label=works+on&message=my+machine&color=green)

My Nextflow Pipelines

## Workflows

## Fastq

  * [Fastqc](workflows/fastqc/) run fastqc on a set of fastq files

### Mapping

  * [Map Fastqs](workflows/mapping/map.fastqs/) Map Fastqs on a reference genome 
  * [Remap bam](workflows/mapping/remap.bwa/) Remap BAM/CRAM on another reference 

### Calling

  * [GATK direct](workflows/gatk/gatk4direct/) Calling bams without gvcfs.
  * [GATK gvcfs](workflows/gatk/gatk4gvcfs/) calling bams with gvcfs.
  * [VQSR](workflows/vqsr/) variant recalibration.

### Structural variants

  * [Indexcov](workflows/indexcov/) Finds CNVs using [indexcov](https://github.com/brentp/goleft/tree/master/indexcov#indexcov).
  * [Delly2](workflows/delly2/) Finds CNVs/SV using [delly](https://github.com/dellytools/delly).
  * [Manta](workflows/manta/) Finds CNVs/SV using using [Manta](https://github.com/Illumina/manta).
  * [Smoove](workflows/smoove/) Finds CNVs/SV using [Smoove](https://github.com/brentp/smoove)
  * [Retrocopies GTF/VCF](workflows/retrocopy/) Finds Retrocopies using a SV VCF and a GTF file.
  * [Plot CNV from VCF](workflows/cnvplotter) Plot CNV as SVG+HTML.

### Relatedness

  * [Somalier BAMS](workflows/somalier/somalier.bams/) Somalier on bams.

### Misc

  * [VCF stats](workflows/vcfstats/vcfstats01/README.md). Apply `bcftools stats` on a set of vcf files
  * [Mosdepth](workflows/mosdepth/). Apply `mosdepth` on a set of BAMs.
  * [FlagStats](workflows/flagstats). Apply 'samtools flagstats' on a set of BAMS.

# Works on my machine

> Why don't you use CONDA ?

I'd love to, but I've got some problems with NFS.

> Why don't you use the `module` directive ?

I will

## About the environment

 * the test for genome is hg19 or hg38 is  here: [./modules/utils/functions.nf](./modules/utils/functions.nf). For now I'm just looking at the fasta filename
 * module `jvarkit` loads java and declares an env variable `${JVARKIT_DIST}` which is the path to the dist directory for the jvarkit jar
 * module `picard` load java and declares an env variable `${PICARD_JAR}` which is the path to picard.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France.

