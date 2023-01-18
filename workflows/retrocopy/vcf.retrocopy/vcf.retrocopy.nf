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
params.vcf = ""
params.gtf = ""
params.publishDir = ""
params.prefix = ""
params.help = false
params.bams = "NO_FILE"

include { VCF_RETROCOPY_01 } from '../../../subworkflows/retrocopy/vcf.retrocopy.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete} from '../../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About

find retrocopies from a SV VCF using a GTF file

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf (file) the indexed vcf file
  * --gtf (file) tabix indexed gtf file
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume vcf.retrocopy.nf  \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcf input.vcf.gz \
	--gtf data.gtf.gz
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
	ch1_ch = VCF_RETROCOPY_01(params, params.reference, params.vcf, params.gtf, file(params.bams))
	html = VERSION_TO_HTML(params,ch1_ch.version)	
	PUBLISH(params,ch1_ch.vcf,html.html,ch1_ch.version,ch1_ch.pdf)
	}


process PUBLISH {
executor "local"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
	val(vcf)
	val(html)
	val(xml)
	val(pdf)
output:
	path("${params.prefix}retrocopy.vcf.gz")
	path("${params.prefix}retrocopy.html")
	path("${params.prefix}retrocopy.xml")
script:
"""
${moduleLoad("bcftools")}

bcftools view -O z -o "${params.prefix}retrocopy.vcf.gz" "${vcf}"
ln -s "${html}" ./${params.prefix}retrocopy.html
ln -s "${xml}" ./${params.prefix}retrocopy.xml
"""
}

runOnComplete(workflow);
