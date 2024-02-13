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

/** path to indexed fasta reference */

include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SEX_GUESS_02} from '../../subworkflows/sex/sex.guess.02.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {	
	ch1 = SEX_GUESS_02(file(params.bams))
	html = VERSION_TO_HTML(ch1.version)
	}



runOnComplete(workflow);

