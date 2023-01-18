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
params.reference = ""
params.fastqs = -1
params.bams = ""
params.help = false
params.publishDir = ""
params.prefix = ""
params.with_bqsr = true

include {MAP_BWA_01} from '../../../subworkflows/mapping/map.bwa.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete} from '../../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

map fastqs on a reference genome

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --fastqs (file) one file containing the paths to the BAM/CRAM. Header: 'sample(tab)R1(tab)R2' [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume map.fastqs.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/referenceX.fasta \\
	--fastqs /path/to/input.tsv
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
	fastqs_ch = Channel.fromPath(params.fastqs).splitCsv(header:true,sep:'\t')
	remap_ch = MAP_BWA_01(params,params.reference,fastqs_ch)
	html = VERSION_TO_HTML(params,remap_ch.version)	
	}


runOnComplete(workflow);

