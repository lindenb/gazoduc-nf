m4_define(`MACRO_HASH', m4_changecom()#m4_changecom(#))m4_dnl
m4_define(`MACRO_PRE', m4_changequote([,])[m4_changequote([,])```[$1]m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`MACRO_CODE', m4_changequote([,])[m4_changequote([,])`[$1]`m4_changequote(`,')]m4_changequote(`,'))m4_dnl
m4_define(`MACRO_DOLLAR',`$'$1)m4_dnl
m4_define(`MACRO_H1',m4_changecom()#m4_changecom(#) $1
)
m4_define(`MACRO_H2',m4_changecom()##m4_changecom(#) $1
)
m4_define(`MACRO_H3',m4_changecom()###m4_changecom(#) $1
)m4_dnl
m4_define(`MACRO_QUOTE',

> $1

)m4_dnl
m4_define(`MACRO_MAKE',`MACRO_PRE(make)
$1
MACRO_PRE')m4_dnl
m4_define(`MACRO_BASH',`MACRO_PRE(bash)
$1
MACRO_PRE')m4_dnl
m4_define(`MACRO_HOW_TO_INSTALL',`

MACRO_H2(System Requirement)

You ll need an install of nextflow. On the Bird cluster, it s available via MACRO_CODE(module)

As I m lazy, the current workflows use local modules, local softwares that may change etc... so it could break your instance of running workflows espcially if it is an old version.
In the future every should move to  MACRO_CODE(docker) MACRO_CODE(singularity) etc...


MACRO_BASH(module purge
module load nextflow)

MACRO_H2(Installation)

If it was not already done clone the repo.

MACRO_BASH(git clone https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf.git)

or if you re using github

MACRO_BASH(git clone https://github.com/lindenb/gazoduc-nf.git)

If the git repository was already installed update the code if needed

MACRO_BASH(cd gazoduc-nf
git pull origin master)

It s also a good practice that you clone the [repo](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) if you need to modify it.

Also if your lab notebook is git-based (of course it is!) you should keep MACRO_CODE(gazoduc-nf) as a git [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules).


On **BiRDCluster** you can also define the following env variables.
MACRO_BASH(export CAPSULE_CACHE_DIR=/LAB-DATA/BiRD/users/`$'{USER}/.nextflow/capsule)



')m4_dnl
m4_dnl
m4_dnl
m4_dnl
m4_define(`MACRO_AWK_PARAMS',m4_changequote([[,]])[[m4_syscmd(awk -vPREFIX="$2" '/M4_END_PARAM/ {exit 0;} {if(match($`0',/^[ \t]*\/[\*]+(.*)[\*]+\/[ \t]*`$'/,a)) {comment=a[1]} else if(match($`0',/^[ \t]*([a-zA-Z][a-zA-Z_0-9]*)[ \t]*=[ \t]*(.*)/,a)>0) { printf("| --%s%s | %s | %s | %s |\n",PREFIX,a[1],a[2],comment, FILENAME); comment="";}  else {comment="";}}' "$1")]]m4_changequote(`,'))m4_dnl
m4_dnl
m4_dnl
m4_dnl
m4_define(`MACRO_OPTIONS',

Parameters can be specified on the command line by prefixing the parameter name with a double dash character e.g. --input.

| Parameters | Default Value | Description | Defined in  |
| --- | --- | --- |
$1

)m4_dnl
m4_dnl
m4_dnl
m4_dnl
m4_define(`MACRO_FOOTER',`

MACRO_H2(Ouput)


Output files are usually available under MACRO_CODE({params.publishDir}/results).

MACRO_H2(On Error)

TODO

MACRO_H2(Workflow)

![workflow.svg](workflow.svg)

MACRO_H2(Author)

 + Pierre Lindenbaum PhD Institut du Thorax. Nantes. France.

')m4_dnl
m4_define(`MACRO_EXECUTE',`

The workflow can be executed using the following command.

MACRO_BASH(module purge

module load nextflow

nextflow -c "../../gazoduc-nf/confs/${HOSTNAME}.cfg" \
        -c /path/to/MACRO_MAIN_CONFIG \
	run \
	-resume \
	-work-dir "/SCRATCH-BIRD/users/${USER}/WORKDIR/" \
	MACRO_MAIN_NF \
	$1
)

where 

 - MACRO_CODE(./../gazoduc-nf/confs/${HOSTNAME}.cfg) is a config allowing to run the SGE job manager on the BirdCluster or SLURM on CCIPL.
 - MACRO_CODE(/path/to/MACRO_MAIN_CONFIG): is the config file containing all the parameters
 - MACRO_CODE(-resume) tell nextflow NOT to restart everything from scratch
 - MACRO_CODE(/SCRATCH-BIRD/users/${USER}/WORKDIR/) is the directory where nextflow will produce the results.
 - MACRO_CODE(MACRO_MAIN_NF) is the main nextflow script




')m4_dnl
m4_define(`MACRO_COMMON_PARAMS',        --genomeId hs37d5 \
        --prefix "20230906.projectName.hs37d5." \
        --publishDir "/SCRATCH-BIRD/users/${USER}/work/projectName/PUBLISH")m4_dnl
m4_dnl
m4_define(`MACRO_PEDIGREE',`MACRO_H2(Pedigree)

A **pedigree** is a tab delimited file without header with the following columns
 
 * family
 * individual (should match the sample name in a VCF  if any)
 * father-id or 0
 * mother-id or 0
 * sex use  MACRO_CODE(male)  or  MACRO_CODE(female) or  MACRO_CODE(0)
 * phenotype# use MACRO_CODE(case) or MACRO_CODE(control) or MACRO_CODE(0)

')m4_dnl
