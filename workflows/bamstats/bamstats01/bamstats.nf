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

include {isBlank;runOnComplete} from '../../../modules/utils/functions.nf'
include {SAMTOOLS_STATS_01} from '../../../subworkflows/samtools/samtools.stats.01.nf'

params.reference=""
params.references="NO_FILE"
params.bams="NO_FILE"
params.bed="NO_FILE"
params.prefix=""
params.publishDir=""
params.help=false


def helpMessage() {
  log.info"""
## About

samtools stats on mulitple files

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAMs [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume bamstats.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also


"""
}


if( params.help ) {
    helpMessage()
    exit 0
}



workflow {
	ch1 = SAMTOOLS_STATS_01(params, params.reference, file(params.references), file(params.bams))
	//PUBLISH(ch1.zip)
	}

process PUBLISH {
tag "${zip.name}"
input:
	path(zip)
when:
	!isBlank(params.getOrDefault("publishDir",""))
script:
"""
mkdir -p "${params.publishDir}"
cp -v "${zip}" "${params.publishDir}/"
"""
}


runOnComplete(workflow);
