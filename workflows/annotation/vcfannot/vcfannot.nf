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
params.prefix = ""
params.genes =""
params.dbsnp =""
params.pedigree =""
params.references=""
params.vcf=""


params.disableFeatures = "tommo"
params.select = ""
params.soacn = "SO:0001574,SO:0001575,SO:0001818"
params.norm  =  false
params.maxalt =  -1
params.description = "no description"
params.gnomadGenome  = ""
params.gnomadExome  = "" 
params.gnomadPop = "AF_nfe"
params.gnomadAF = 0.001
params.vcfid = ""
params.mergevcf = "NO_FILE"
params.joinvcf = ""
params.joinsuffix = "_EXTERNAL"
params.extraBcfTools = ""
params.extension = ".vcf.gz"
params.hardfilters = ""
params.lowDP = 10
params.lowGQ = 70
params.captures = ""
params.gnomadSVAF = 0.1
params.goTerms = "GO:0006936,GO:0022857,GO:0060047,GO:0003179,GO:0008016"
params.lowMQ = 30
params.mqRankSum = 5
params.readPosRankSum = 5
params.sor = 5
params.polyx = 10
/** limit to this bed */
params.bed = "NO_FILE"

include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete} from '../../../modules/utils/functions.nf'
include {ANNOTATE_VCF_02} from '../../../subworkflows/annotation/annotation.vcf.02.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'

def helpMessage() {
  log.info"""
## About


## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcf path/to/vcf
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
	ann_ch = ANNOTATE_VCF_02(params, params.reference, file(params.vcf), file(params.bed))

	html = VERSION_TO_HTML(params, ann_ch.version)
	
	//publish_ch = Channel.empty().mix(html.html).mix(version_ch).mix(concat_ch.vcf).mix(concat_ch.index)
	
	//PUBLISH(publish_ch.collect())
	}


process PUBLISH {
tag "N=${L.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(L)
output:
        path("*.bcf"),optional:true
        path("*.bcf.csi"),optional:true
        path("*.vcf.gz"),optional:true
        path("*.vcf.gz.tbi"),optional:true
        path("*.xml")
        path("*.html")
when:
        !params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
for F in ${L.join(" ")}
do
        ln -s "\${F}" ./
done
"""
}

runOnComplete(workflow);
