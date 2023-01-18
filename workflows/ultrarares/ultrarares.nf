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
params.bams = "NO_FILE"
params.vcf = "NO_FILE"
/* run ultrare for each record in the BED file OR make the bed by default */
params.bed = "NO_FILE"
params.gnomad_max_af=0.01
params.gnomad_population="AF_nfe"
params.n_bams_per_hc_call=5
params.help = false
params.publishDir = ""
params.prefix = ""
params.extraBcftoolsView1=""
params.extraBcftoolsView2=""

include {ULTRA_RARES_01} from '../../subworkflows/ultrarares/ultrarares.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

search for ultra rares variants

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf (file) indexed vcf file. default: ""
  * --bams (file) controls bams
  * --bed (file) how to split the vcf
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcf /path/to/in.vcf.gz \\
	--bams /path/to/bams.list \\
	--bed /path/input.bed
```

## Workflow

![workflow](./workflow.svg)
  

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	ch1 = ULTRA_RARES_01(params, params.reference, params.vcf,file(params.bams), file(params.bed))
	html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.vcf,ch1.index,ch1.pdf)
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(version)
	val(html)
	val(vcf)
output:
	path("*.bcf")
	path("*.csi")
	path("*.html")
	path("*.xml")
	path("*.pdf")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${vcf} ./
ln -s ${csi} ./
ln -s ${html} ./
ln -s ${version} ./
ln -s ${pdf} ./
"""
}

runOnComplete(workflow);

