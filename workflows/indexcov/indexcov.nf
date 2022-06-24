/*

Copyright (c) 2022 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
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
/** files prefix */
params.prefix = ""
/** goleft version */
params.goleft_version = "v0.2.4"

include {INDEXCOV} from '../../subworkflows/indexcov/indexcov.nf'
include {assertFileExists} from '../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

Detects CNV using go-left indexcov

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --mapq (int)  min mapping quality . If it's lower than 0 (this is the default) just use the bam index as is. Otherwise, rebuild the bai
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --goleft_version (string) default: "${params.goleft_version}"

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume indexcov.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--mapq 30
```

## Workflow

![workflow](./workflow.svg)
  
## See also

 * indexcov: https://github.com/brentp/goleft/tree/master/indexcov
 * https://twitter.com/yokofakun/status/1527419449669734426

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	assertFileExists(params.reference,"--reference must be defined")
	assertFileExists(params.bams,"--bams must be defined")
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

