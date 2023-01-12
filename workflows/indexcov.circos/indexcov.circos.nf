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

gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()

gazoduc.build("bed","NO_FILE").argName("bed").desc("The BED file generated by indexcov").required().put()

include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'
include {INDEXCOV_TO_CIRCOS_01} from '../../subworkflows/indexcov/indexcov.circos.01.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

if(params.help) {
	gazoduc.usage().name("indexcov2circos").desc("Plot circos from indexcov output").print();
	exit 0
	}
else
	{
	gazoduc.validate();
	}



workflow {
	ch1 = INDEXCOV_TO_CIRCOS_01(params,params.reference,file(params.bed))
	html = VERSION_TO_HTML(params, ch1.version)
	to_zip =  Channel.empty().mix(ch1.version).mix(html.html).mix(ch1.png).mix(ch1.svg)
	SIMPLE_ZIP_01(params, to_zip.collect())
	}

runOnComplete(workflow)
