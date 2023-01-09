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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults();

/** path to indexed fasta reference */
params.reference = ""
/** one file containing the paths to the BAM/CRAM  */
params.bams = ""
params.with_merge_manta_vcf = false
params.manta_cpus = 16

include {runOnComplete} from '../../modules/utils/functions.nf'
include {MANTA_SINGLE_SV01} from '../../subworkflows/manta/manta.single.01.nf' 
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'


if( params.help ) {
    gazoduc.usage().
	name("manta").
	description("call SV/CNV with manta for each bam").
	print()
    exit 0
} else {
	gazoduc.validate();
	}


workflow {
	manta_ch = MANTA_SINGLE_SV01(params, params.reference, params.bams)


	html = VERSION_TO_HTML(params,manta_ch.version)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(manta_ch.version).
			mix(manta_ch.zip).
			mix(manta_ch.merge_vcf).
			mix(manta_ch.merge_vcf_index).
			mix(manta_ch.manta_files).
			mix(html.html)


	PUBLISH(to_publish.collect())
	}

runOnComplete(workflow)

process PUBLISH {
tag "N=${files.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(files)
output:
	path("*.bcf")
	path("*.bcf.csi")
	path("*.xml")
	path("*.html")
	path("*.zip")
	path("*.list")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = params.prefix?:""
"""
for F in ${files.join(" ")}
do
	ln -s "\${F}" ./
done
"""
}

