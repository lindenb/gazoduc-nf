# IGV Report

generate static HTML report from VCF + BAM using https://github.com/igvteam/igv-reports.

## Author(s)

  * Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.

## Options

### Main options

    * --bams <path to file> . File containing the paths to the BAM/CRAMS files. One path per line. [REQUIRED]. [NO_FILE].
    * --num_cases <integer> . number of samples carrying a ALT allele to choose from the VCF. [5].
    * --num_controls <integer> . number of samples carrying *NO* ALT allele to choose from the VCF. [5].
    * --vcf <path to file> . file containing the variant. [REQUIRED]. [NO_FILE].


### Help

    * --help <true|false> . Display help for this workflow and exit. [true].

### Input

    * --reference <path to fasta> . Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict' ( or picard CreateSequenceDictionary ) . [REQUIRED]. [false].

### Output

    * --prefix <string> . set a suffix for the files generated for this workflow. [].
    * --publishDir <directory> . set a base directory where final output files should be written.. [].


## Issues

report issues at https://github.com/lindenb/gazoduc-nf/issues


