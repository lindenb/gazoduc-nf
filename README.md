GAZODUC
=======

![Last commit](https://img.shields.io/github/last-commit/lindenb/gazoduc-nf.png)
![works on my machine](https://img.shields.io/static/v1?label=works+on&message=my+machine&color=green)
[![pipeline status](https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf/badges/master/pipeline.svg)](https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf/-/commits/master)

My [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) [Nextflow](https://www.nextflow.io) Pipelines used at the [Insitut du Thorax/Nantes/France](https://umr1087.univ-nantes.fr/).

Works on my machine


## Installation

**Glicid Users** : you can try to install everything using conda/micromamba but I had a poor experience with this solution


## Install java


Install a recent version of openjdk (e.g: jdk-23.0.1) from https://jdk.java.net/23/ , and for example download and expand it into ${HOME}/packages/jdk-23.0.1 ( or look under https://jdk.java.net/archive/ ).


add JDK to your ${PATH}:
```
export PATH=${HOME}/packages/jdk-23.0.1/bin:${PATH}
```

set the variable JAVA_HOME
```
export JAVA_HOME=${HOME}/packages/jdk-23.0.1
```
## Install conda/micromamba

**Glicid Users** : see https://doc.glicid.fr/GLiCID-PUBLIC/alpha/logiciels.html#_installer_micromamba_sur_glicid

## Install nextflow

```
$ mkdir -p ${HOME}/packages
cd ${HOME}/packages
curl -s https://get.nextflow.io | bash
# make it executable
chmod +x nextflow
```

create the directories where nextflow will download its libraries:

```
mkdir -p "/scratch/nautilus/users/${USER}/.nextflow/capsule"
```

**Glicid Users** :  Prepare a CONDA directory for nextflow:

```
mkdir -p /micromamba/${USER}/envs/nextflow-envs
```

**Glicid Users** :   Open ${HOME}/.bash_profile and add the following lines. Adjust the variable according to your installation

```
export NAUTILUS_SCRATCH=/scratch/nautilus/users/${USER}
# disable nextflow logging because I want one line per job instead of one line per process
export NXF_ANSI_LOG=false
# where to put conda/micromamba stuff
export NXF_CONDA_CACHEDIR=/micromamba/${USER}/envs/nextflow-envs
# if you want to use nf-core, don't connect to the web to check the version
export NFCORE_NO_VERSION_CHECK=1
# nextflow wants to put file here
export CAPSULE_CACHE_DIR=${NAUTILUS_SCRATCH}/.nextflow/capsule
# override NF home in the scratch instead of your home
export NXF_HOME=${NAUTILUS_SCRATCH}/.nextflow
# where to install singularity/apptainer stuff
export NXF_SINGULARITY_CACHEDIR=${NAUTILUS_SCRATCH}/.singularity/cache
export NXF_APPTAINER_CACHEDIR=${NXF_SINGULARITY_CACHEDIR}
export NXF_SINGULARITY_LIBRARYDIR=${NAUTILUS_SCRATCH}/.singularity/lib
export NXF_APPTAINER_LIBRARYDIR=${NXF_SINGULARITY_LIBRARYDIR}
# set java home
export JAVA_HOME=${HOME}/packages/jdk-23.0.1
export PATH=${JAVA_HOME}/bin:${HOME}/packages/nextflow:${PATH}
# the following line was found to be required if you want to use micromamba with nextflow
export PATH=${HOME}/.local/bin:${PATH}
```


Reload `bash`, Test your installation of nextflow:

```
$ nextflow info
```

## nextflow Profiles for **GLICID** users


https://nextflow.io/docs/latest/config.html#config-profiles

> Configuration files can define one or more profiles. A profile is a set of configuration settings that can be selected during pipeline execution using the -profile command line option.



 * `micromamba` : tell nextflow to use micromamba
 * `GRCh38` : fills required information for the GRCh38 reference
 * `hs38me` : fills required information for the hs38me reference
 * `hs37d5` : fills required information for the hs37d5 reference
 * `waves` : tell nextflow to use the `waves` cluster  .`export SLURM_CLUSTERS=waves` should be set.
 * `nautilus` : tell nextflow to use the `nautilus` cluster .`export SLURM_CLUSTERS=nautilus` should be set.


Example:

```
nextflow run -profile 'GRCh38,nautilus` main.nf
```

## Install GAZODUC

```
git clone https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf.git
```


## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France.
