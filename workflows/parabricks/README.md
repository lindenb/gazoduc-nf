# Parabricks

Parabricks Pipeline

## Input / Ouput

General Parameters

| Parameter | Description | Type | Default | Required | Pattern |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA REFERENCE file. | `string` |  | True | ^\S+\.(fasta|fa|fna)?(\.gz)?$ |
| `prefix` | prefix for the save files. A good prefix would be 'YYYY-MM-DD.project_name.'  <details><summary>Help</summary><small>Ouput file prefix</small></details>| `string` |  |  | ^[a-zA-z0-9_\.]+\.$ |
| `hts_type` | Type of analysis (WES or WGS). If WES, --capture (bed) is expected. <details><summary>Help</summary><small>Type of HTS analysis</small></details>| `string` |  |  | ^(WES|WGS)$ |
| `bed` | Path to BED file that was used as the capture. Required for --hts_type=WES. <details><summary>Help</summary><small>Capture BED</small></details>| `string` |  |  |  |
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
