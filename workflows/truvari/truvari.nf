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

gazoduc.build("vcfs","NO_FILE").
	desc("file containing the path to the VCF files. One per line.").
	existingFile().
	required().
	put()




include {TRUVARI_01} from '../../subworkflows/truvari/truvari.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("truvari").
	desc("merge vcfs using truvari https://github.com/ACEnglish/truvari").
	print();
    exit 0
} else {
   gazoduc.validate();
}


workflow {
	ch1 = TRUVARI_01(params,params.reference,file(params.vcfs))
	html = VERSION_TO_HTML(params,ch1.version)

	pub_ch = Channel.empty().mix(ch1.vcf).mix(ch1.version).mix(html.html)
	SIMPLE_PUBLISH_01(params, pub_ch.collect())
	}

runOnComplete(workflow)
