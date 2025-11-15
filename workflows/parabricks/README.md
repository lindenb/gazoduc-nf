# Parabricks

Parabricks Pipeline

## Input / Ouput

General Parameters

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA REFERENCE file. | `string` |  | True |
| `prefix` | prefix for the save files. A good prefix would be 'YYYY-MM-DD.project_name.'  <details><summary>Help</summary><small>Ouput file prefix</small></details>| `string` |  |  |
| `hts_type` | Type of analysis (WES or WGS). If WES, --capture (bed) is expected. <details><summary>Help</summary><small>Type of HTS analysis</small></details>| `string` |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. <details><summary>Help</summary><small>output directory</small></details>| `string` |  |  |
| `help` |  | `boolean` |  |  |

## Mapping FASTQs

Parameters use to map the reads on a reference genome.

| Parameter | Description | Type | Default | Required |
|-----------|-----------|-----------|-----------|-----------|
| `samplesheet` | Path to samplesheet with the appropriate suffix (csv or tsv or json). <details><summary>Help</summary><small>path to fastq samplesheet.</small></details>| `string` |  |  |
| `orad_directory` | ora-d directory containing the sofwtare used to decrompress .ora files. If not defined, the toolbox is downloaded. (But on glicid it's blocked by the proxy) <details><summary>Help</summary><small>path to ora-d directory used to decrompress .ora files.</small></details>| `string` |  |  |
| `mapper` | wich mapper to use. bwa,pb,parabricks,dragmap | `string` |  | True |
| `known_indels_vcf` | Path to  VCF indexed with tabix used as source of known indels for BQSR. | `string` |  |  |
| `with_merge_fastqs` | If true, merge the reads for the same sample. Always true for mapper=parabricks | `boolean` |  |  |
| `with_fastp` | If true, preprocess the Reads with FASTP before mapping the reads | `boolean` |  |  |
| `with_fastqc` | Run Fastqc. If --with_fastp==true, anothe fastqc is run after trimming. | `boolean` |  |  |
| `with_cram` | Save the alignments as CRAM (always true for parabricks). | `boolean` |  |  |
| `bam2fastq_method` | when the input is a BAM/CRAM to has be converted back to fastq prior to mapping, what is the method to use ?  | `string` |  |  |
