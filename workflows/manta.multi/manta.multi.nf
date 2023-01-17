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
/** one file containing the paths to the BAM/CRAM  */
params.bams = ""
/** display help */
params.help = false
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.prefix = ""
params.bed = "NO_FILE"

include {runOnComplete} from '../../modules/utils/functions.nf'
include {MANTA_MULTI_SV01} from '../../subworkflows/manta/manta.multi.01.nf' 
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

def helpMessage() {
  log.info"""
## About

Detects CNV/SV using manta, joint calling

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM. [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

  * https://github.com/Illumina/manta

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	manta_ch = MANTA_MULTI_SV01(params, params.reference, file(params.bams), file(params.bed))

	PUBLISH(manta_ch.zip)
	}

process PUBLISH {
executor "local"
tag "${zip.name}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(zip)
output:
	path("*.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s "${zip.toRealPath()}" ./
"""
}

runOnComplete(workflow)
