# Parabricks pipeline parameters

Parabricks Pipeline

## Input / Ouput

General Parameters

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA REFERENCE file. | `string` |  | True | ^\S+\.(fasta|fa|fna)?(\.gz)?$ |
| `prefix` | prefix for the save files. A good prefix would be 'YYYY-MM-DD.project_name.'  <details><summary>Help</summary><small>Ouput file prefix</small></details>| `string` |  |  | ^[a-zA-z0-9_\.]+\.$ |
| `hts_type` | Type of analysis (WES or WGS). If WES, --capture (bed) is expected. <details><summary>Help</summary><small>Type of HTS analysis</small></details>| `string` |  |  | ^(WES|WGS)$ |
| `bed` | Path to BED file that was used as the capture. Required for --hts_type=WES. For WGS a good idea is to make a BED excluding the hard-to-sequence regions (e.g: https://github.com/Boyle-Lab/Blacklist/  ) <details><summary>Help</summary><small>Capture BED</small></details>| `string` |  |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. <details><summary>Help</summary><small>output directory</small></details>| `string` |  |  |  |
| `help` |  | `boolean` |  |  |  |

## Mapping FASTQs

Parameters use to map the reads on a reference genome.

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` | Path to samplesheet with the appropriate suffix (csv or tsv or json).<br/>For **ora** files or **fastq.gz** from illumina/dragen, the 'sample' column is optional as it can be extracted from the fastq name <details><summary>Help</summary><small>path to fastq samplesheet.</small></details>| `string` |  |  | ^\S+\.(csv|tsv|json)$ |
| `orad_directory` | ora-d directory containing the sofwtare used to decrompress .ora files. If not defined, the toolbox is downloaded. (But on glicid it's blocked by the proxy) <details><summary>Help</summary><small>path to ora-d directory used to decrompress .ora files.</small></details>| `string` |  |  |  |
| `mapper` | wich mapper to use. one of 'bwa,pb,parabricks,dragmap' | `string` |  | True | ^(bwa|pb|parabricks|dragmap)$ |
| `known_indels_vcf` | Path to  VCF indexed with tabix used as source of known indels for BQSR. | `string` |  |  | ^\S+\.vcf\.gz$ |
| `with_merge_fastqs` | If true, merge the reads for the same sample. Always true for mapper=parabricks | `boolean` |  |  |  |
| `with_fastp` | If true, preprocess the Reads with FASTP before mapping the reads | `boolean` |  |  |  |
| `with_cram` | Save the alignments as CRAM (always true for parabricks). | `boolean` |  |  |  |
| `bam2fastq_method` | when the input is a BAM/CRAM to has be converted back to fastq prior to mapping, what is the method to use ?  | `string` |  |  | ^(samtools|picard|gatk|pb|parabricks)$ |

## QC for BAMs

Parameters used for BAM QC

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_fastqc` | Run Fastqc. If --with_fastp==true, anothe fastqc is run after trimming. | `boolean` |  |  |  |
| `with_collect_metrics` | Run Picard Collect Metrics (which can be quite slow for WGS). | `boolean` |  |  |  |
| `mosdepth_use_median` | use 'median' instead of 'average' to calculate the depth with mosdepth. | `boolean` |  |  |  |
| `mosdepth_min_depth` | Treshold of total mosdepth coverage. Bam below this treshold withh be affected by '--on_low_bam_quality'. | `number` |  |  |  |
| `on_low_bam_quality` | What to do with BAMs with a depth under the treshold ? abort: stop the workflow. warning : just print a message. skip: skip the bams with a low quality. | `string` |  |  | ^(abort|warning|skip)$ |

## SOMALIER

Parameters used for SOMALIER

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_somalier` | Run **somalier** on BAMs. | `boolean` |  |  |  |
| `somalier_vcf_sites` | somalier VCF sites. if ignored, VCF is downloaded . | `string` |  |  |  |

## MANTA

Manta: Structural variant and indel caller for mapped sequencing data  (https://github.com/Illumina/manta)

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_manta` | Run manta. | `boolean` |  |  |  |
| `with_truvari` | Apply truvari to merge the VCFs (https://github.com/ACEnglish/truvari) | `boolean` |  |  |  |

## SMOOVE

Smoove: 'Structural variant calling and genotyping with existing tools, but, smoothly'.  https://github.com/brentp/smoove

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_smoove` | Run Smoove. | `boolean` |  |  |  |

## DELLY2

 DELLY2: Structural variant discovery by integrated paired-end and split-read analysis . https://github.com/dellytools/delly

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_delly` | Run DELLY2. | `boolean` |  |  |  |

## INDEXCOV

 Goleft/indexcov: Quickly estimate coverage from a whole-genome bam or cram index.  https://github.com/brentp/goleft/tree/master/indexcov/

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_indexcov` | Run Goleft/indexcov. | `boolean` |  |  |  |
| `indexcov_batchsize` | group BAMs by batch of 'x' bams. Useful when manupulating thousands of bams | `number` |  |  |  |

## CNVnator

CNVnator: a tool for CNV discovery and genotyping from depth-of-coverage by mapped reads  : https://github.com/abyzovlab/CNVnator

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_cnvnator` | Run CNVnator | `boolean` |  |  |  |
| `cnvnator_bin_size` | CNVnator 'bin' size | `number` |  |  |  |

## SNV Calling

General parameters for the SNV calling

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `makewindows_args` | Arguments for 'betools makewindows' Before calling the SNVs, the genome or the capture will be divided into smaller parts. For exome this parameter is quite useless unless the window size is lower than the fragment size. See also --bed_cluster_args2  | `string` |  |  |  |
| `bed_cluster_args2` | Arguments for 'jvarkit bed cluster' group the fragments generated using --makewindows_args into bed file of size 'x'. For example for exomes, you can imagine grouping the fragments by group of 'x' bp. Used by graphtyper. | `string` |  |  |  |
| `exclude_encode_blacklist` | Exclude encode balcklist regions from the calling regions | `boolean` |  |  |  |

## Graphtyper

Population-scale genotyping using pangenome graphs :  https://github.com/DecodeGenetics/graphtyper

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_graphtyper` | Run graphtyper | `boolean` |  |  |  |

## Freebayes

Freebayes :  Bayesian haplotype-based genetic polymorphism discovery and genotyping. https://github.com/freebayes/freebayes 

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_freebayes` | Run freebayes | `boolean` |  |  |  |
| `freebayes_args` | extra parametes for freebayes | `string` |  |  |  |

## Bcftools

Bcftools  call variants with bcftools . https://github.com/samtools/bcftools 

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_bcftools_call` | Run bcftools call | `boolean` |  |  |  |

## GATK HaplotypeCaller

GATK Haplotypecaller in GVCF mode. 

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_haplotype_caller` | Run GATK Haplotypcaller | `boolean` |  |  |  |
| `hc_gvcf_merge_method` | Method to merge and genotype GVCF files ? 'combinegvcfs' or 'glnexus'  | `string` |  |  | ^(combinegvcfs|glnexus)$ |

## Parabricks DeepVariant

Deep Variant :  DeepVariant is an analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data.  https://github.com/google/deepvariant

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_pb_deepvariant` | Run Parabricks DeepVariant | `boolean` |  |  |  |

## Parabricks HaplotypeCaller

Run the parabricks version of haplotypecaller.

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_pb_haplotypecaller` | Run Parabricks haplotypecaller | `boolean` |  |  |  |

## De Novo assembly of unmapped reads

Assemble unmapped reads with spades

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `with_unmapped` | try to assemble unmapped reads with spades | `boolean` |  |  |  |
| `unmapped_fast_mode` | fast mode: use the last reads in the bam with unmapped ('*') chromosome. Slow mode: scan the whole bam, extracting discordant read in pairs where one read is unmapped. | `boolean` |  |  |  |
