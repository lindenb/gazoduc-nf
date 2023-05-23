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


gazoduc.build("bams", "NO_FILE").
	desc("File containing the paths to the indexed BAM/CRAM files.").
	existingFile().
	required().
	put()

gazoduc.build("gtf", "NO_FILE").
	desc("gtf file indexed with tabix").
	put()

gazoduc.build("bed", "NO_FILE").
	desc("SV intervals as a BED file").
	existingFile().
	required().
	put()

gazoduc.build("mapq","0").
	desc("min mapping quality").
	setInt().
	put()

gazoduc.build("median",true).
	desc("normalize on median").
	setBoolean().
	put()

gazoduc.build("extend_bed",3.0).
	desc("extend original bed by this factor").
	setDouble().
	put()

gazoduc.build("max_sv_length",-1).
	desc("max abs(SV_LENGTH) or <0 to ignore").
	setInt().
	put()


params.references="NO_FILE"

include {PLOT_COVERAGE_01} from '../../subworkflows/plotdepth/plot.coverage.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


if( params.help ) {
    gazoduc.usage().
	name("plot coverage").
	desc("Plot coverage for a set of bams").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch1 = PLOT_COVERAGE_01(params,params.reference, file(params.references),file(params.bams), file(params.bed), file(params.gtf))
	//html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.zip)
	}

runOnComplete(workflow);


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


