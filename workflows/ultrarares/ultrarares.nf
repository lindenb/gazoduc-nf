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


def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
	required().
	existingFile().
        put()

gazoduc.make("vcf","NO_FILE").
        description("initial vcf file or use --bams").
        put()

gazoduc.make("bed","NO_FILE").
        description("run ultrare for each record in the BED file OR make the bed by default").
        put()

gazoduc.make("mapq",-1).
        description("mapping quality or -1").
        setInt().
        put()

gazoduc.make("gnomad_max_af",0.001).
        description("gnomad max AF").
        setDouble().
        put()

gazoduc.make("n_bams_per_hc_call",30).
        description("number of BAMS per gatk HC call").
        setInt().
        put()


gazoduc.make("vcf2interval_distance","50mb").
        description("split VCF per region of 'x' size").
        put()


params.gnomad_population="AF_nfe"
params.extraBcftoolsView1=""
params.extraBcftoolsView2=""

include {ULTRA_RARES_01} from '../../subworkflows/ultrarares/ultrarares.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


if( params.help ) {
    gazoduc.usage().
        name("Ultrarares").
        desc("ultrares").
        print();
    exit 0
} else {
   gazoduc.validate();
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

