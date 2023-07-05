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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults()

gazoduc.make("vcf","NO_FILE").
	description("path to an indexed VCF or BCF file").
	required().
	existingFile().
	put()

gazoduc.make("bams","NO_FILE").
	description("file containing the path to multiple bam files").
	required().
	existingFile().
	put()

gazoduc.make("mapq",0).
	description("min mapping quality").
	setInt().
	put()

gazoduc.make("excludeids","NO_FILE").
	description("file containing IDS to exclude").
	put()



include {CNV_PLOTTER_01} from '../../subworkflows/cnvplotter/cnv.plotter.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

if( params.help ) {
    gazoduc.usage().
		name("cnvplotter").
		description("plot CNVs as SVG.").
		print();
    exit 0
    }
else
	{
	gazoduc.validate()
	}


workflow {
	ch1 = CNV_PLOTTER_01(params.genomeId, params.vcf , params.bams, file(params.excludeids))
	html = VERSION_TO_HTML(ch1.version)
	PUBLISH(ch1.version,html.html,ch1.zip)
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(version)
	val(html)
	val(zip)
output:
	path("*.zip")
	path("*.html")
	path("*.xml")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${version} ./
ln -s ${html} ./
ln -s ${zip} ./
"""
}

runOnComplete(workflow);
