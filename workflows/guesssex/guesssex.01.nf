/*

Copyright (c) 2024 Pierre Lindenbaum

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
params.prefix = ""
params.publishDir = ""
params.bams = "NO_FILE"
params.help=false

include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SEX_GUESS_02} from '../../subworkflows/sex/sex.guess.02.nf'

def helpMessage() {
  log.info"""
## About

guess samples' sex from BAM files.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) path to multiple bams files [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams
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
	ch1 = SEX_GUESS_02(params, params.reference, file(params.references), file(params.bams))
	html = VERSION_TO_HTML(params,ch1.version)
	ch2 =  Channel.empty().mix(ch1.version).mix(html.html).mix(ch1.pdf).mix(ch1.output)
	PUBLISH(params,ch2.collect())
	}


process PUBLISH {
executor "local"
publishDir "${meta.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
	val(L)
when:
	meta.containsKey("publishDir")
script:
"""
cp -v ${L.join(" ")} "${meta.publishDir}"
"""
}


runOnComplete(workflow);

