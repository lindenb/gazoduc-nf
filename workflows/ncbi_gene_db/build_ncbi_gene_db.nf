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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults();


include {runOnComplete} from '../../modules/utils/functions.nf'
include {BUILD_NCBI_GENE_DB} from '../../subworkflows/ncbi/build_ncbi_gene_db.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("build ncbi db").
	description("build table including ncbi genename/geneid/genedescription").
	print()
    exit 0
} else {
	gazoduc.validate();
	}


workflow {
    ch1 = BUILD_NCBI_GENE_DB(params)
    html = VERSION_TO_HTML(params,ch1.version)

    SIMPLE_PUBLISH_01(params, Channel.empty().mix(html.html).mix(ch1.version).mix(ch1.output).collect())
    }

runOnComplete(workflow)

