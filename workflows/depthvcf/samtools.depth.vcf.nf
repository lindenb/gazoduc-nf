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
params.references = "NO_FILE"
params.mapq = 1
params.bams = ""
params.bed = "NO_FILE"
params.help = false
params.publishDir = ""
params.prefix = ""

include {SAMTOOLS_DEPTH_VCF_01} from '../../subworkflows/samtools/samtools.depth.vcf.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


def helpMessage() {
  log.info"""
## About

apply samtools depth to a set of bams.

## Author

Pierre Lindenbaum

## Options
  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --references (file) path to known references.
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --bed (file) optional BED file. default: ""
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--bed /path/to/in.bed
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
	ch1 = SAMTOOLS_DEPTH_VCF_01(params,params.reference,file(params.references), file(params.bams), file(params.bed))
	html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.output)
	}


process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(version)
	path(html)
	path(summary)
output:
	path("${params.prefix}samtools.depth.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
zip -j "${params.prefix}samtools.depth.zip" ${version} ${html} ${summary}
"""
}

runOnComplete(workflow);

