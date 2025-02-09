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

include {PLOT_COVERAGE_01} from '../../subworkflows/plotdepth'
include {dumpParams; runOnComplete} from '../../modules/utils/functions.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	reference = Channel.of(file(params.fasta), file(params.fai), file(params.dict)).collect()
	ch1 = Channel.fromPath(params.bams).
		splitText().
		map{it.trim()}.
		map{[file(it),file(it+(it.endsWith(".bam")?".bai":".crai"))]}

	ch1 = PLOT_COVERAGE_01(reference, ch1 , file(params.bed))
	}

runOnComplete(workflow);

