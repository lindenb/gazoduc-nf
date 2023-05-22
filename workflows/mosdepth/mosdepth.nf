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

gazoduc.make("mapq",1).
	description("mapq min mapping quality . If it's <=0, just use the bam index as is. Otherwise OR IF IT's A CRAM, rebuild the bai").
	setInteger().
	argName("MAPQ").
	put()

gazoduc.make("bams","NO_FILE").
	description("File containing the paths to the BAM/CRAMS files. One path per line").
	required().
	existingFile().
	put()

gazoduc.make("mosdepth_extra","").
	description("Extra parameters for mosdepth").
	put()

gazoduc.make("bed","NO_FILE").
	description("limit to that BED file").
	put()



include {MOSDEPTH_BAMS_01} from '../../subworkflows/mosdepth/mosdepth.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

if( params.help ) {
    gazoduc.usage().
		name("mosdepth").
		description("apply mosdepth to a set of bams").
		print();
    exit 0
    } else {
	gazoduc.validate();
	}


workflow {
	ch1 = MOSDEPTH_BAMS_01(params,params.reference,params.bams, file(params.bed))
	html = VERSION_TO_HTML(params,ch1.version)
	PUBLISH(ch1.version,html.html,ch1.summary,ch1.multiqc,ch1.pdf.collect())
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(version)
	path(html)
	path(summary)
	path(multiqc)
	path(pdfs)
output:
	path("${params.prefix}mosdepth.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
zip -j "${params.prefix}mosdepth.zip" ${version} ${html} ${summary} ${multiqc} ${pdfs.join(" ")}
"""
}

runOnComplete(workflow);

