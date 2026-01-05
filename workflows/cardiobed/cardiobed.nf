/*

Copyright (c) 2026 Pierre Lindenbaum

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
params.gtf = ""
params.goTerms = ""
params.excludeGoTerms = ""
params.help = false
params.publishDir = ""
params.prefix = ""

include {runOnComplete} from '../../modules/utils/functions.nf'
include {CARDIOBED_01} from '../../subworkflows/go/cardiobed.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'

def helpMessage() {
  log.info"""
## About

create Bed for cardiac genes

## Author

Pierre Lindenbaum PhD

## Options

  * --reference (fasta) Indexed Fasta Reference file [REQUIRED]
  * --gtf (file) GTF file
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--gtf /path/to/gtf
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
	ch1 = CARDIOBED_01(params, params.reference, file(params.gtf))
	html = VERSION_TO_HTML(params, ch1.version)
	pub_ch = Channel.empty().mix(html.html).mix(ch1.version).mix(ch1.bed).mix(ch1.index).mix(ch1.terms)
	SIMPLE_PUBLISH_01(params, pub_ch.collect())
	}

runOnComplete(workflow)

