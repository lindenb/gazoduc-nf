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
/** mapq min mapping quality . If it's <=0, just use the bam index as is. Otherwise, rebuild the bai */
params.mapq = -1
/** one file containing the paths to the BAM/CRAM */
params.bams = ""
/** display help */
params.help = false
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.references=""
params.prefix = ""
params.dbsnp=""
params.pedigree=""
params.beds=""
params.nbams=20

include {GATK4_HAPCALLER_DIRECT_01} from '../../../subworkflows/gatk/gatk4.hapcaller.direct.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad} from '../../../modules/utils/functions.nf'

def helpMessage() {
  log.info"""
## About


## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --beds (file) call in the following bed files. One path to bed per line.
  * --mapq (int)  min mapping quality . If it's lower than 0 (this is the default) just use the bam index as is. Otherwise, rebuild the bai
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume gatk4direct.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--beds /path/to/beds.list \\
	--mapq 30
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
	publish_ch = Channel.empty()

	gatk_ch= GATK4_HAPCALLER_DIRECT_01(
			params,
			params.reference,
			params.bams,
			file(params.beds)
			)
	html_ch = VERSION_TO_HTML(params,gatk_ch.version)

	publish_ch = publish_ch.mix(gatk_ch.version)
	publish_ch = publish_ch.mix(html_ch.html)

	PUBLISH(gatk_ch.vcf, gatk_ch.index, publish_ch.collect())
	}


process PUBLISH {
tag "${L.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(vcf)
	path(index)
	val(L)
output:
	path("*.bcf"),optional:true
	path("*.csi"),optional:true
	path("*.xml"),optional:true
	path("*.html"),optional:true
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
cp "${vcf}" "${params.prefix?:""}genotyped.bcf"
cp "${index}" "${params.prefix?:""}genotyped.bcf.csi"

for F in ${L.join(" ")}
do
        ln -s "\${F}" ./
done
"""
}

runOnComplete(workflow)
