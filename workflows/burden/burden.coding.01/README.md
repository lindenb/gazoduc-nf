# burden

burden for coding regions

## Author(s)

  * Pierre Lindenbaum PhD. Institut du Thorax. U1087. 44400 Nantes. France.

## Options

### Main options

    * --bed <value> . limit analysis to that BED file. [NO_FILE].
    * --bed_cluster_method <value> . when clustering exons with jvarkit/bedcluster, how do we cluster genes.. [--size 5mb].
    * --genes_slop <integer> . extend exons by 'x' bases. [50].
    * --pedigree <path to file> . Jvarkit formatted pedigree. Tab delimited, no header, FAM/ID/FATHER/MOTHER/SEX/PHENOTYPE. [REQUIRED]. [NO_FILE].
    * --vcf <value> . Path to a VCF file or a file with the .list' suffix containing the full path to several VCFs file. The VCF must be indexed with 'bcftools index' (an associated .tbi/.csi must be present) . [REQUIRED]. [NO_FILE].
    * --with_pihat <true|false> . Run a Pihat before burden. [false].


### Gnomad

    * --gnomad_exome_hg19 <value> . gnomad exome path for hg19. [/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/release-181127/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.vcf.gz].
    * --gnomad_genome_hg19 <value> . gnomad genome path for hg19. [/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/release-181127/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.gz].
    * --gnomad_genome_hg38 <value> . gnomad genome path for hg38. [/LAB-DATA/BiRD/resources/species/human/broadinstitute.org/gnomad/3.0/gnomad.genomes.r3.0.sites.vcf.gz].

### Help

    * --help <true|false> . Display help for this workflow and exit. [true].

### Input

    * --reference <path to fasta> . Path to the reference genome as FASTA. The file must be indexed with 'samtools faidx' and 'samtools dict' ( or picard CreateSequenceDictionary ) . [REQUIRED]. [false].

### Output

    * --prefix <string> . set a suffix for the files generated for this workflow. [].
    * --publishDir <directory> . set a base directory where final output files should be written.. [].

### pihat

    * --pihat_MAF <double> . pihat MAF. [0.1].
    * --pihat_f_missing <double> . pihat fraction of missing genotypes. [0.05].
    * --pihat_filters <value> . pihat filter. [ --apply-filters '.,PASS' ].
    * --pihat_max <double> . max pihat . [0.05].
    * --pihat_max_DP <integer> . max Depth. [300].
    * --pihat_min_DP <integer> . min Depth. [10].
    * --pihat_min_GQ <integer> . pihat min GQ. [20].

### rvtest

    * --rvtest_arguments <value> . Parameters for rvtest. This should include the tests to be performed. [REQUIRED]. [--burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat --vt price,analytic --kernel 'skat[nPerm=1000],kbac,skato'].

### wgselect

    * --minRatioSingleton <double> . remove variant if HET singleton has AD ratio out of x< AD/ratio < (1.0-x). [0.2].
    * --wgselect_MQ <double> . INFO/RMSMappingQuality (MQ) It is meant to include the standard deviation of the mapping qualities. Including the standard deviation allows us to include the variation in the dataset. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants. Dicard variant with MQ<x. Ignore if value is <= 0.. [10].
    * --wgselect_MQRankSum <double> . This is the u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping qualities of the reads supporting the reference allele and the alternate allele. A positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele; a negative value indicates the mapping qualities of the reference allele are higher than those supporting the alternate allele. A value close to zero is best and indicates little difference between the mapping qualities See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants.Dicard variant with MQRankSum<x Ignore if value is >0.. [-10].
    * --wgselect_ReadPosRankSum <double> . INFO/ReadPosRankSum compares whether the positions of the reference and alternate alleles are different within the reads. See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants. Dicard variant with -x<ReadPosRankSum<x. Ignore if value is 0.. [4].
    * --wgselect_annot_method <value> . how to annotate ? 'vep' or 'snpeff'. [vep].
    * --wgselect_cadd_phred <double> . Discard variants having CADD phred treshold < 'x'. Ignore if 'x' < 0.0 or --wgselect_cadd_tabix is not defined.. [-1.0].
    * --wgselect_cadd_tabix <value> . path to tabix-indexed fill containing CADD data. (like ... CADD/v1.6/whole_genome_SNVs.tsv.gz ). [].
    * --wgselect_exclude_soacn <value> . Accesion number for consequences to keep after the functional annotation. If not empty, variant with the following SO will be removed after the standard selection on --wgselect_soacn .Multiple are comma separated. [].
    * --wgselect_f_missing <double> . Discard variants where the proportion of NO_CALL (./.) genotypes is grater than is value. [0.05].
    * --wgselect_fisherh <double> . remove variant if fisher test per variant is lower than 'x'. Disable if <0.. [0.05].
    * --wgselect_gnomadAF <0.0 - 1.0> . Allele frequency(AF) treshold in GNOMAD VCF. Discard variants with INFO/AF> value. [0.01].
    * --wgselect_gnomadPop <population> . Allele frequency(AF) to consider in GNOMAD VCF. [AF_nfe].
    * --wgselect_hwe <double> . remove variants with HW test. Ask Floriane :-P . Disable if <0.. [1E-15].
    * --wgselect_inverse_so <true|false> . Inverse SO. [false].
    * --wgselect_lowGQ <integer> . ALL genotypes carrying a ALT must have a Genotype Quality GQ >= x. Ignore if x <=0 . [70].
    * --wgselect_mapability_hg19_bigwig <path to file> . Mapability for hg19. bigwig file with the number of time a region is repeated in the genome. 1=uniq . discard if !=1. [/LAB-DATA/BiRD/resources/species/human/ucsc/hg19/encodeDCC/wgEncodeDukeMapabilityUniqueness35bp.bigWig].
    * --wgselect_mapability_hg38_bigwig <path to file> . Mapability for hg38. bigwig file with the number of time a region is repeated in the genome. 1=uniq . discard if !=1. [/LAB-DATA/BiRD/resources/species/human/ucsc/hg38/hoffmanMappability/k24.Umap.MultiTrackMappability.bw].
    * --wgselect_maxDP <double> . remove variant mean called genotype depth is tool high.. [300].
    * --wgselect_max_alleles_count <integer> . select variants on min/max number of alleles (diallelic is 2). [3].
    * --wgselect_maxmaf <double> . remove variant if internal MAF is too high. Disable if < 0. [0.05].
    * --wgselect_minDP <double> . remove variant mean called genotype depth is tool low.. [10].
    * --wgselect_minGQsingleton <integer> . remove variant if singleton has bad GQ < x. [99].
    * --wgselect_minRatioSingleton <double> . remove variant if HET singleton has AD ratio out of x< AD/ratio < (1.0-x). [0.2].
    * --wgselect_polyx <size> . remove variant near a poly-x with size > value. [10].
    * --wgselect_soacn <value> . Accesion number for consequences to keep after the functional annotation. Multiple are comma separated. [SO:0001574,SO:0001575,SO:0001818].
    * --wgselect_with_contrast <true|false> . Apply bcftools contrast on VCF. [true].
    * --wgselect_with_count <true|false> . Count variants at each step of wgselect. [true].
    * --wgselect_with_encode_exclude <true|false> . remove variant overlapping excluded ENCODe regions. [true].
    * --wgselect_with_homvar <true|false> . remove variant on autosome if no HET and found at least one HOM_VAR. [true].
    * --wgselect_with_lcr <true|false> . remove variant overlapping low complexity regions. [true].
    * --wgselect_with_rmsk <true|false> . remove variant overlapping repeat masker UCSC regions. [true].
    * --wgselect_with_simple_repeats <true|false> . remove variant overlapping simple repeats regions. [true].


## Issues

report issues at https://github.com/lindenb/gazoduc-nf/issues


