# trackhub

Build a UCSC track hub form a set of SV/CNV VCFs files.

## Author(s)

  * Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.

## Options

### Main options

    * --description <value> . track description. [REQUIRED]. [].
    * --max_sv_len <integer> . max SV length. [10000000].
    * --min_sv_len <integer> . min SV length. [0].
    * --vcfs <path to file> . File with the .list' suffix containing the full path to several VCFs file. The VCF must be indexed with 'bcftools index' (an associated .tbi/.csi must be present) . [REQUIRED]. [NO_FILE].
    * --with_bnd <true|false> . include variants with INFO/SVTYPE=BND. [true].


### Genomes

    * --genomeId <value> . The main genome used. This is the genome id in the XML file (see option --genomesFile). [REQUIRED]. [false].
    * --genomesFile <path to file> . Path to a XML file describing all the available genomes on your server. See doc. [REQUIRED]. [false].

### Help

    * --help <true|false> . Display help for this workflow and exit. [true].

### Output

    * --prefix <string> . set a suffix for the files generated for this workflow. [].
    * --publishDir <directory> . set a base directory where final output files should be written.. [].


## Issues

report issues at https://github.com/lindenb/gazoduc-nf/issues


