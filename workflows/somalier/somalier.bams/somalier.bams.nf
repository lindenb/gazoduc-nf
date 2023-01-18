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
params.bams=""
params.pedigree=""
params.help = false
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.prefix = ""

include {SOMALIER_BAMS_01} from  '../../../subworkflows/somalier/somalier.bams.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete} from '../../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About


## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --mapq (int)  min mapping quality . If it's lower than 0 (this is the default) just use the bam index as is. Otherwise, rebuild the bai
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume somalier.bams.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

* https://github.com/brentp/somalier

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}

workflow {

	somalier_ch = SOMALIER_BAMS_01(params,params.reference,
		Channel.fromPath(params.bams),
		file(params.pedigree.isEmpty()?"NO_FILE":params.pedigree)
		)

        html = VERSION_TO_HTML(params,somalier_ch.version)
	PUBLISH(somalier_ch.zip, somalier_ch.version, html.html)
	}

process PUBLISH {
tag "${zip.name}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(zip)
	val(version)
	val(html)
output:
	path("*.zip")
	path("*.html")
	path("*.xml")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
for F in ${zip} ${version} ${html}
do
	ln -s "\${F}" ./
done
"""
}

runOnComplete(workflow)
