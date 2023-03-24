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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
        required().
        existingFile().
        put()

gazoduc.make("mapq",-1).
        description("mapping quality or -1").
        setInt().
        put()


gazoduc.make("nbams",20).
        description("number of bams per HC").
        setInt().
        put()

gazoduc.make("references","").
        description("other references").
        put()

gazoduc.make("dbsnp","").
        description("path to dbsnp").
        put()

gazoduc.make("pedigree","").
        description("path to pedigree").
        put()

gazoduc.make("beds","").
        description("path to a list of bed files").
        put()


include {GATK4_HAPCALLER_DIRECT_01} from '../../../subworkflows/gatk/gatk4.hapcaller.direct.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad} from '../../../modules/utils/functions.nf'

if( params.help ) {
    gazoduc.usage().
        name("GraphTyper CNV").
        desc("Genotype CNV/SV using graphtyper").
        print();
    exit 0
} else {
   gazoduc.validate();
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
