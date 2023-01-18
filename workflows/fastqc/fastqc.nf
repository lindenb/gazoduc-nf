/*

Copyright (c) 2023 Pierre Lindenbaum

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
params.fastqs = ""
params.help = false
params.publishDir = ""
params.prefix = ""

include {FASTQC_01} from '../../subworkflows/fastqc/fastqc.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

run FASTC for a list of bams

## Author

Pierre Lindenbaum

## Options

  * --fastqs (file) list of fastqs one per line [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume fastqc.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--fastqs paths.list
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
	ch1 = FASTQC_01(params,Channel.fromPath(params.fastqs))
	PUBLISH(ch1.zip)
	}

process PUBLISH {
executor "local"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(zip)
output:
	path("*.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${zip} ./
"""
}

runOnComplete(workflow);

