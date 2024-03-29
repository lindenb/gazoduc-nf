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
params.references= "NO_FILE"
params.bams = ""
params.bed = "NO_FILE"
params.help = false
params.publishDir = ""
params.prefix = ""

include {BAM_FLAGSTATS_01} from '../../subworkflows/samtools/samtools.flagstats.01.nf' 
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

apply samtools flagstats for a set of bams

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --references 
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume mosdepth.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--references /path/to/references.txt \\
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
	ch1 = BAM_FLAGSTATS_01(params,params.reference,file(params.references),file(params.bams))
	html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.summary,ch1.pdf.collect())
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(version)
	path(html)
	path(summary)
	path(pdfs)
output:
	path("${params.prefix}mosdepth.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
zip -j "${params.prefix}mosdepth.zip" ${version} ${html} ${summary} ${pdfs.join(" ")}
"""
}

runOnComplete(workflow);

