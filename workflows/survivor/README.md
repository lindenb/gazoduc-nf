# survivor

merge vcfs using survivor https://github.com/fritzsedlazeck/SURVIVOR

## Author(s)

  * Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.

## Options

### Main options

    * --vcfs <path to file> . file containing the path to the VCF files. One per line.. [REQUIRED]. [NO_FILE].


### Help

    * --help <true|false> . Display help for this workflow and exit. [true].

### Input

    * --reference <path to fasta> . Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict' ( or picard CreateSequenceDictionary ) . [REQUIRED]. [false].

### Output

    * --prefix <string> . set a suffix for the files generated for this workflow. [].
    * --publishDir <directory> . set a base directory where final output files should be written.. [].

### Survivor

    * --survivor_merge_params <value> . Parameters for survivor (max distance between breakpoints ; Minimum number of supporting caller; Take the type into account;  Take the strands of SVs into accoun; Disabled; Minimum size of SVs to be taken into account).. [500 1 1 1 0 100].


## Issues

report issues at https://github.com/lindenb/gazoduc-nf/issues


