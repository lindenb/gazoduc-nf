GAZODUC
=======

![Last commit](https://img.shields.io/github/last-commit/lindenb/gazoduc-nf.png)
![works on my machine](https://img.shields.io/static/v1?label=works+on&message=my+machine&color=green)

My Nextflow Pipelines

## Workflows

  * [Indexcov](workflows/indexcov/README.md) Finds CNVs using [indexcov](https://github.com/brentp/goleft/tree/master/indexcov#indexcov).
  * [Delly2](workflows/delly2/README.md) Finds CNVs/SV using [delly](https://github.com/dellytools/delly).
  * [Manta](workflows/manta/README.md) Finds CNVs/SV using using [Manta](https://github.com/Illumina/manta).
  * [Smoove](workflows/smoove/README.md) Finds CNVs/SV using [Smoove](https://github.com/brentp/smoove)

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

