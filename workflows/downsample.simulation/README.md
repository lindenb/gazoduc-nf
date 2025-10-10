# Motivation

downsamples BAM/CRAM files, generates downsamples bam files, statistics about GT concordance

# Installation

See [Installation](../../README.md]

# parameters

## Required

  * `--fasta` full path to FASTA reference ( or use profile GRCh38 on glicid)
  * `--samplesheet` : CSV file. Required columns
    - 'bam': full path to BAM file

## Optional

  * `--depths` : comma separated or wanted file. default is `1,2,5,10,15,20`
  * `--bed`: limit to that BED file.
  * `--prefix` : files prefix
  * `--outdir`:  where to write results.

## Example

```
echo -e 'chr3\t38522656\t38675094' > gene.bed
echo "bam" > samplesheet.csv && find DIR -type f -name "*.bam" >> samplesheet.csv

export SLURM_CLUSTERS=waves && \
nextflow run  -resume -work-dir work /path/to/gazoduc-nf/workflows/downsample.simulation/main.nf run \
        -profile "micromamba,${SLURM_CLUSTERS},GRCh38" \
        --prefix "20251009.downsample." \
        --outdir "results" \
        --depths "5,15" \
        --samplesheet samplesheet.csv \
        --bed gene.bed

```
## Output

output folder contains each bam downsampled, the bam index, the outputs of mosdepth.
multiqc directory,  and pipeline infos from nextflow.

```
/results/
|-- 20251009.downsample.GRCh38.lowpass.pdf
|-- BAMS
|   |-- 01-055-C
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP001.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP001.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP001.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP002.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP002.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP002.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP005.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP005.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP005.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP010.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP010.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP010.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP015.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP015.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP015.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP020.bam
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP020.bam.bai
|   |   |-- 20251009.downsample.GRCh38.01-055-C.DP020.bam.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.gene.checked.vcf.gz
|   |   |-- 20251009.downsample.GRCh38.01-055-C.gene.checked.vcf.gz.md5
|   |   |-- 20251009.downsample.GRCh38.01-055-C.gene.checked.vcf.gz.tbi
|   |   |-- 20251009.downsample.GRCh38.01-055-C.mosdepth.global.dist.txt
|   |   |-- 20251009.downsample.GRCh38.01-055-C.mosdepth.region.dist.txt
|   |   |-- 20251009.downsample.GRCh38.01-055-C.mosdepth.summary.txt
|   |   |-- 20251009.downsample.GRCh38.01-055-C.regions.bed.gz
|   |   `-- 20251009.downsample.GRCh38.01-055-C.regions.bed.gz.csi
(...)
|-- multiqc
|   `-- all
|       `-- 20251009.downsample.GRCh38.multiqc.zip
`-- pipeline_info
    |-- execution_report_2025-10-10_16-55-02.html
    |-- execution_report_2025-10-10_16-55-44.html
    |-- execution_report_2025-10-10_16-57-50.html
```
