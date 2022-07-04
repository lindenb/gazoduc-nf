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

include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {moduleLoad;readContigFile} from '../../../modules/utils/functions.nf'
include {VCF_TO_BED_01} from '../../../modules/jvarkit/jvarkit.vcf2bed.01.nf'
include {ANNOTATE} from '../../../subworkflows/annotation/annotation.vcf.01.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'

def helpMessage() {
  log.info"""
## About


## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --genes (file) gene list, one gene per line.

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume gatk4genes.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--genes /path/to/genes.txt
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
	version_ch = Channel.empty()

	c1_ch = VCF_TO_BED_01(params,params.vcf)
	version_ch = version_ch.mix(c1_ch.version)

	cfg2=readContigFile(params.annotation_config?:"${workflow.projectDir}/../../../confs/annotation.cfg").plus(params)


	annotate_ch = ANNOTATE(cfg2 , params.reference, c1_ch.output.splitCsv(header:true,sep:'\t'))
	version_ch = version_ch.mix(annotate_ch.version)

	
	file_list_ch = COLLECT_TO_FILE_01([:],annotate_ch.bedvcf.
			splitCsv(header:false,sep:'\t').
			map{T->T[1]}.collect())
	
	concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output)
	version_ch = version_ch.mix(concat_ch.version)

	
	version_ch = MERGE_VERSION(params, "annot vcf", "annotation of VCF", version_ch.collect())
	html = VERSION_TO_HTML(params,version_ch.version)
	
	publish_ch = Channel.empty().mix(html.html).mix(version_ch).mix(concat_ch.vcf).mix(concat_ch.index)
	
	PUBLISH(publish_ch.collect())
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

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

