/*

Copyright (c) 2025 Pierre Lindenbaum

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

include {runOnComplete;dumpParams    } from '../../modules/utils/functions.nf'
include {MANTA_GERMLINE_SINGLE_SV01  } from '../../subworkflows/manta/single' 
include {PREPARE_REFERENCE           } from '../../subworkflows/samtools/prepare.ref'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	version_ch = Channel.empty()
	def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name: (params.ucsc_name?:"undefined")
		]
	def fasta = [ hash_ref, file(params.fasta)]
	PREPARE_REFERENCE(hash_ref,fasta)
	version_ch = version_ch.mix(PREPARE_REFERENCE.out.versions)

	ch0 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{assertKeyExistsAndNotEmpty(it,"sample")}
        .map{h->hasKey(h,"id")?h:h.plus(id:h.sample)}


	manta_ch = MANTA_GERMLINE_SINGLE_SV01(
		hash_ref,
		fasta,
		PREPARE_REFERENCE.out.fai,
		PREPARE_REFERENCE.out.dict,
		Channel.fromPath(params.samplesheet),
		file(params.bed)
		)
	}
runOnComplete(workflow)


