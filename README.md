GAZODUC
=======

![Last commit](https://img.shields.io/github/last-commit/lindenb/gazoduc-nf.png)
![works on my machine](https://img.shields.io/static/v1?label=works+on&message=my+machine&color=green)

My [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) [Nextflow](https://www.nextflow.io) Pipelines used at the [Insitut du Thorax/Nantes/France](https://umr1087.univ-nantes.fr/).

Warning: Some NF files use local configurations like `module`, `conda`, full paths to resources. Some workflows are missing an example, tests, etc...

Warning: some NF use a custom java library under `lib/`. The NF workflow should be invoked with `nextflow run -lib /path/to/lib ... `

## Workflows

## Fastq

  * [Fastqc](workflows/fastqc/) run fastqc on a set of fastq files

### Mapping

  * [Map Fastqs](workflows/mapping/map.fastqs/) Map Fastqs on a reference genome 
  * [Remap bam](workflows/mapping/remap.bwa/) Remap BAM/CRAM on another reference 
  * [Fastmap](workflows/mapping/fastmap/) Fast map Fastq (split fastq, bwa-mem, samtools merge)

### Calling

  * [GATK direct](workflows/gatk/gatk4direct/) Calling bams without gvcfs.
  * [GATK gvcfs](workflows/gatk/gatk4gvcfs/) calling bams with gvcfs.
  * [VQSR](workflows/vqsr/) variant recalibration.
  * [Graphtyper](workflows/graphtyper/graphtyper.genotype/) calling with graphtyper
  * [Deepvariant](workflows/deepvariant) calling with DeepVariant
  * [bcftools](workflows/bcftools/calling) calling with bcftools/mpileup+call
  * [ultrarares](workflows/ultrarares). Call variants from a VCF for a list a BAM.

### Annotation

  * [VCF fastannot](annotation/fastannot/) quick VCF annotation
  * [GO / BED](workflows/go4bed) generate BED file from GFF3+GO+GOA
  * [CardioBED](workflows/cardiobed) generate BED file from GFF3+GO+GOA for Cardiac-related genes

### Structural variants / CNV

#### Calling

  * [Indexcov](workflows/indexcov/) Finds CNVs using [indexcov](https://github.com/brentp/goleft/tree/master/indexcov#indexcov).
  * [Delly2](workflows/delly2/) Finds CNVs/SV using [delly](https://github.com/dellytools/delly).
  * [Manta](workflows/manta/) Finds CNVs/SV using using [Manta](https://github.com/Illumina/manta). One VCF per sample.
  * [Manta-Multi](workflows/manta.multi/) Finds CNVs/SV using using [Manta](https://github.com/Illumina/manta). Join calling.
  * [Smoove](workflows/smoove/) Finds CNVs/SV using [Smoove](https://github.com/brentp/smoove)
  * [Retrocopies GTF/VCF](workflows/retrocopy/vcf.retrocopy) Finds Retrocopies using a SV VCF and a GTF file.
  * [Duphold](workflows/duphold) Apply Duphold to a set of VCFs+BAMs.
  * [etching](workflows/etching.germline) Call germline SV using etching. Doesn't work for now.
  * [graphtyper CNV](workflows/graphtyper/graphtyper.cnv) Genotype SV variants using graphtype.

#### Merging

  * [Sruvivor](workflows/survivor) Merge SV using survivor.
  * [truvari](workflows/truvari) Merge SV using truvari.

#### Visualization

  * [Plot CNV from VCF](workflows/cnvplotter) Plot CNV as SVG+HTML.
  * [Plot CNV using bed+ samtools depth](workflows/plotcoverage) Plot CNV as PDF.
  * [indexcov+circos](workflows/indexcov.circos) Plot indexcov output as circos
  * [wgs coverage](workflows/wgscovplot) WGS coverage genome wide. Plot as SVG.
  * [UCSC track hub](workflows/trackhub) Generate Track Hub for ucsc from CNV VCFs

### RNAseq

  * [Sashimi/Encode](workflows/sashimi.encode/) sashimi plots from Encode data.

### Viruses

  * [Gridss/Virus BreakEnd]workflows/gridss/virusbreakend/).  GRIDSS viral integration detection tool (the vcf was empty for the few BAMs I tested)

### Relatedness

  * [Somalier BAMS](workflows/somalier/somalier.bams/) Somalier on bams.
  * [Somalier VCFs](workflows/somalier/somalier.vcfs/) Somalier on VCFs.

### Burden

  * [Burden Coding](workflows/burden/burden.coding.01/) Burden on coding regions.
  * [Burden First Intron](workflows/burden/burden.first.intron/) Burden on first intron.
  * [Burden UTR](workflows/burden/burden.utr/) Burden on UTR5', UTR3', UTR5+3'.
  * [Burden uORF](workflows/burden/burden.uorf.vep/) Burden on micro-ORF (uORF) in UTR. Using VEP plugin.
  * [Burden Sliding windows](burden/burden.sliding.window/) Burden with sliding windows.
  * [Burden Pairs](burden/burden.pairs) Burden per pairs of genes.
  * [WGSelect](workflows/wgselect/basic) Basic "wgselect".
  * [Optimize](workflows/burden/optimize.rvtests) Find peak of p-value by applying sliding window on a gene.

### VNTR

  * [ExpansionHunter](workflows/expansionhunter/) Apply Illumina/ExpansionHunter to a set of bams.
  * [ExpansionHunter DeNovo](workflows/expansionhunter.denovo)

### Misc

  * [VCF stats](workflows/vcfstats/vcfstats01/README.md). Apply `bcftools stats` on a set of vcf files
  * [Samtools stats](workflows/bamstats/bamstats01/README.md). Apply `samtools stats` on a set of bam files
  * [Mosdepth](workflows/mosdepth/). Apply `mosdepth` on a set of BAMs.
  * [Samtools depth](workflows/samtoolsdepth/). Apply `samtools depth` on a set of BAMs.
  * [Samtools depth VCF](workflows/depthvcf/). Apply `samtools depth` on a set of BAMs. ouput is a vcf file.
  * [FlagStats](workflows/flagstats/). Apply 'samtools flagstats' on a set of BAMS.
  * [SpliceAI](workflows/spliceai). Apply SpliceAI to a VCF
  * [Sex](workflows/guesssex) . Guess sex from chrX/chrY in BAM file
  * [QC remapping](workflows/comparemapping). QC remapping by genotype concordance + liftover
  * [contimination](workflows/contamination). Contaminations in BAMs/FASTQs.
  * [IGV Reports](workflows/igvreport). Generate IGV static HTML reports from VCF+BAMS.

### Biostars

Solutions for biostars.org.

  + [Biostars9550377](workflows/biostars/biostars9550377.nf) speeding up bcftools view.



# Works on my machine

> Why don't you use CONDA ?

I'd love to, but I've got some problems with NFS.

> Why don't you use/write NF-CORE ?

same: I've got some problems with NFS.

> Why don't you use the `module` directive ?

I will... I will

## About the environment

 * the test for genome is hg19 or hg38 is  here: [./modules/utils/functions.nf](./modules/utils/functions.nf). For now I'm just looking at the fasta filename
 * module `jvarkit` loads java and declares an env variable `${JVARKIT_DIST}` which is the path to the dist directory for the jvarkit jar
 * module `picard` load java and declares an env variable `${PICARD_JAR}` which is the path to picard.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France.
