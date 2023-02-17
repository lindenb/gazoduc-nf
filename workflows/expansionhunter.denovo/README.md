# expansion hunter de novo

expansion hunter de novo

## Author(s)

  * Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.

## Options

### Main options

    * --cases <path to file> . file containing the path to multiple bam files for the cases. [REQUIRED]. [NO_FILE].
    * --controls <path to file> . file containing the path to multiple bam files for the controls. [REQUIRED]. [NO_FILE].


### Help

    * --help <true|false> . Display help for this workflow and exit. [true].

### Input

    * --reference <path to fasta> . Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict' ( or picard CreateSequenceDictionary ) . [REQUIRED]. [false].

### Output

    * --prefix <string> . set a suffix for the files generated for this workflow. [].
    * --publishDir <directory> . set a base directory where final output files should be written.. [].


## Issues

report issues at https://github.com/lindenb/gazoduc-nf/issues


