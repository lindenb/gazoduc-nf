


# burdencoding


## About


Burden test for rare variants in coding regions.



## System Requirement


You ll need an install of nextflow. On the Bird cluster, it s available via `module`

As I m lazy, the current workflows use local modules, local softwares that may change etc... so it could break your instance of running workflows espcially if it is an old version.
In the future every should move to  `docker` `singularity` etc...


```bash
module purge
module load nextflow
```

## Installation


If it was not already done clone the repo.

```bash
git clone https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf.git
```

or if you re using github

```bash
git clone https://github.com/lindenb/gazoduc-nf.git
```

If the git repository was already installed update the code if needed

```bash
cd gazoduc-nf
git pull origin master
```

It s also a good practice that you clone the [repo](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) if you need to modify it.

Also if your lab notebook is git-based (of course it is!) you should keep `gazoduc-nf` as a git [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules).


On **BiRDCluster** you can also define the following env variables.
```bash
export CAPSULE_CACHE_DIR=/LAB-DATA/BiRD/users/${USER}/.nextflow/capsule
```





## Pedigree


A **pedigree** is a tab delimited file without header with the following columns
 
 * family
 * individual (should match the sample name in a VCF  if any)
 * father-id or 0
 * mother-id or 0
 * sex use  `male`  or  `female` or  `0`
 * phenotype# use MACRO_CODE(case) or MACRO_CODE(control) or MACRO_CODE(0)



## Parameters


| Field | Description | Value | Source |
| `--prefix` | files will be generated with this prefix .eg: "20230101.hello." |  | `../../../confs/by_workflow/../default.params.cfg` |
| `--help` | should you print the help and exit ? | false | `../../../confs/by_workflow/../default.params.cfg` |
| `--publishDir` | base directory where results will be written |  | `../../../confs/by_workflow/../default.params.cfg` |
| `--genomeId` | genome identifier in gazoduc-nf/confs/params.genomes. e.g: "hs37d5" |  | `../../../confs/by_workflow/../genomeId.params.cfg` |
| `--gatk.mapq` |  | 20 | `../../../confs/by_workflow/../by_subworkflow/../by_app/gatk.cfg` |
| `--gatk.haplotypecaller.mapq` |  | 10 | `../../../confs/by_workflow/../by_subworkflow/../by_app/gatk.cfg` |
| `--gatk.haplotypecaller.maxAlternateAlleles` |  | 6 | `../../../confs/by_workflow/../by_subworkflow/../by_app/gatk.cfg` |
| `--gatk.collectWgsMetrics.args` |  |  MINIMUM_MAPPING_QUALITY=${params.gatk.mapq} MINIMUM_BASE_QUALITY=20 COVERAGE_CAP=250 LOCUS_ACCUMULATION_CAP=100000 STOP_AFTER=100000 INCLUDE_BQ_HISTOGRAM=false | `../../../confs/by_workflow/../by_subworkflow/../by_app/gatk.cfg` |
| `--gatk.intervalListToBed.args` |  |  --SCORE 500 --SORT true | `../../../confs/by_workflow/../by_subworkflow/../by_app/gatk.cfg` |
| `--wgselect.distance` | when splitting vcf into parts, make interval of that distance | 10mb | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.min_distance` | when splitting vcf into parts, don't leave variant if distance lower that this value | 100 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.bcftools_options` | option for first bcftools, e.g: --apply-filters '.,PASS' |  | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.gatk_hardfiltering_percentile` | apply gatk hard filtering ignore if < 0 | 0.001 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.f_missing` | fraction of missing allele | 0.05 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_rmsk` | filter out variant overlapping repeat masker data in ucsc | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_encode_exclude` | filter out variants in encode blacklisted | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_lcr` | filter out variants in low complexity region | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_simple_repeats` | filter out variant overlapping  ucsc  simple repeat | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.max_alleles_count` | max alleles per variant | 3 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.polyx` | max polyx for jvarkit/vcfployx | 10 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_kinship` | use kinship : not sure it is still used | false | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.fisherh` | fisher horizontal : remove variant if fisher test per variant is lower than 'x'. Disable if <0. | 0.05 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.soacn` | keep so consequences | SO:0001629,SO:0001818 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.gnomadPop` | GNOMAD population | AF_nfe | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.gnomadAF` | gnomad max allele frequency | 0.01 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.minGQsingleton` | remove variant if singleton has bad GQ < x | 99 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.lowGQ` | ALL genotypes carrying a ALT must have a Genotype Quality GQ >= x. Ignore if x <=0 | 70 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.annot_method` | how to annotate ? 'vep' or 'snpeff' | snpeff | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.minDP` | remove variant mean called genotype depth is tool low | 10 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.maxDP` | remove variant mean called genotype depth is tool high | 300 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_count` | Count variants at each step of wgselect | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_homvar` | remove variant on autosome if no HET and found at least one HOM_VAR | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.maxmaf` | remove variant if internal MAF is too high. Disable if < 0 | 0.1 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.hwe` | remove variants with HW test. Ask Floriane :-P . Disable if <0. | 0.000000000000001 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.with_contrast` | Apply bcftools contrast on VCF | true | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.inverse_so` | inverse output of vcffilterso | false | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.minRatioSingleton` | emove variant if HET singleton has AD ratio out of x< AD/ratio < (1.0-x) | 0.2 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.cadd_phred` | TODO | -1.0 | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.gnomadgenomefilterexpr_hg19` | remove gnomad hg19 expression | FILTER~\"GNOMAD_GENOME_BAD_AF\"&#124;&#124; FILTER~\"GNOMAD_GENOME_InbreedingCoeff\"&#124;&#124; FILTER~\"GNOMAD_GENOME_RF\" | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--wgselect.gnomadgenomefilterexpr_hg38` | remove gnomad hg38 expression | FILTER~\"GNOMAD_GENOME_BAD_AF\"&#124;&#124; FILTER~\"GNOMAD_GENOME_InbreedingCoeff\"&#124;&#124; FILTER~\"GNOMAD_GENOME_AS_VQSR\" | `../../../confs/by_workflow/../by_subworkflow/wgselect.config` |
| `--vcf` | indexed VCF or a list of vcf with the .list suffix | NO_FILE | `../../../confs/by_workflow/wgselect.basic.cfg` |
| `--pedigree` | pedigree | NO_FILE | `../../../confs/by_workflow/wgselect.basic.cfg` |
| `--bed` | limit to that bed | NO_FILE | `../../../confs/by_workflow/wgselect.basic.cfg` |




The workflow can be executed using the following command.

```bash
module purge

module load nextflow

nextflow -c "../../gazoduc-nf/confs/${HOSTNAME}.cfg" \
        -c /path/to/MACRO_MAIN_CONFIG \
	run \
	-resume \
	-work-dir "/SCRATCH-BIRD/users/${USER}/WORKDIR/" \
	gazoduc-nf/burden/burden.coding.01/burden.coding.01.nf \
	--vcf /path/to/vcf \
	--genomeId hs37d5 \
        --prefix "20230906.projectName.hs37d5." \
        --publishDir "/SCRATCH-BIRD/users/${USER}/work/projectName/PUBLISH" \
	--pedigree "/path/to/file.ped" \
	--bed /path/to/input.bed 


```

where 

 - `./../gazoduc-nf/confs/${HOSTNAME}.cfg` is a config allowing to run the SGE job manager on the BirdCluster or SLURM on CCIPL.
 - `/path/to/MACRO_MAIN_CONFIG`: is the config file containing all the parameters
 - `-resume` tell nextflow NOT to restart everything from scratch
 - `/SCRATCH-BIRD/users/${USER}/WORKDIR/` is the directory where nextflow will produce the results.
 - `gazoduc-nf/burden/burden.coding.01/burden.coding.01.nf` is the main nextflow script








## Ouput



Output files are usually available under `{params.publishDir}/results`.

## On Error


TODO

## Workflow


![workflow.svg](workflow.svg)

## Author


 + Pierre Lindenbaum PhD Institut du Thorax. Nantes. France.


