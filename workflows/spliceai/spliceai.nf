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

params.conda = ""


gazoduc.
    make("vcf","NO_FILE").
    description(gazoduc.Gazoduc.DESC_VCF_OR_VCF_LIST).
    required().
    put()

gazoduc.
    make("bed","NO_FILE").
    description("limit analysis to that BED file").
    put()

gazoduc.
    make("pedigree", "NO_FILE").
    description(gazoduc.Gazoduc.DESC_JVARKIT_PEDIGREE).
    put()

include {SPLICEAI_01} from '../../subworkflows/spliceai/spliceai.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


if(params.help) {
        gazoduc.usage().name("spliceai").description("spliceai").print()
        exit 0
        }
else
        {
        gazoduc.validate()
        }


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	ch1 = SPLICEAI_01([:], params.genomeId, Channel.fromPath(params.vcf), file(params.pedigree), file(params.bed))
	html = VERSION_TO_HTML(ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.summary,ch1.pdf.collect())
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(version)
	path(html)
	path(summary)
	path(pdfs)
output:
	path("${params.prefix}mosdepth.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
zip -j "${params.prefix}mosdepth.zip" ${version} ${html} ${summary} ${pdfs.join(" ")}
"""
}

runOnComplete(workflow);

