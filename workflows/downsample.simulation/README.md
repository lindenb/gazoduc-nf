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
        --samplesheet samplesheet.csv \
        --bed ene.bed

```
