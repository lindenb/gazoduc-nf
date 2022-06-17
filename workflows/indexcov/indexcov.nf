nextflow.enable.dsl=2

/** path to indexed fasta reference */
params.reference = ""
/** mapq min mapping quality . If it's <=0, just use the bam index as is. Otherwise, rebuild the bai */
params.mapq = -1
/** one file containing the paths to the BAM/CRAM */
params.bams = ""
/** display help */
params.help = false
/** publish Directory */
params.publishDir = ""

include {INDEXCOV} from '../../subworkflows/indexcov/indexcov.nf'


def helpMessage() {
  log.info"""
## About

Detects CNV using go-left indexcov

## Author

Pierre Lindenbaum Phd

## Options

  * --reference (fasta) indexed fasta reference [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --mapq (int)  min mapping quality . If it's lower than 0 (this is the default) just use the bam index as is. Otherwise, rebuild the bai
  * --publishDir (dir) Save output in this directory

## Usage

```
module load nextflow && nextflow -C ../../confs/cluster.cfg  run -resume indexcov.nf \
	--publishDir output \
	 -with-trace -with-timeline \
	--reference /LAB-DATA/BiRD/resources/species/human/cng.fr/hs37d5/hs37d5_all_chr.fasta \
	--bams /path/to/bams.list \
	--mapq 30
```
  
"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	indexcov_ch = INDEXCOV(params,params.reference,Channel.fromPath(params.bams))
	PUBLISH(indexcov_ch.zip)
	}

process PUBLISH {
tag "${zip.name}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(zip)
output:
	path(zip)
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
echo "publishing ${zip}"
"""
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

