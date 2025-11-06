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

include { runOnComplete;dumpParams    } from '../../modules/utils/functions.nf'
include { assertKeyExistsAndNotEmpty  } from '../../modules/utils/functions.nf'
include { META_TO_BAMS                } from '../../subworkflows/samtools/meta2bams1'
include { MANTA_SINGLE                } from '../../subworkflows/manta/single' 
include { PREPARE_ONE_REFERENCE       } from '../../subworkflows/samtools/prepare.one.ref'
if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	versions = Channel.empty()
	def metadata = [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name: (params.ucsc_name?:"undefined")
		]
	def fasta = [ metadata, file(params.fasta)]
	PREPARE_ONE_REFERENCE(metadata,Channem.of(fasta))
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	if(params.samplesheet.endsWith(".json")) {
		ch0 = Channel.fromPath(params.samplesheet)
        	.splitJSon()
		}
	else
		{
		ch0 = Channel.fromPath(params.samplesheet)
        	.splitCsv(header:true,sep:',')
		}

	META_TO_BAMS(
  			metadata,
  			PREPARE_ONE_REFERENCE.out.fasta,
  			PREPARE_ONE_REFERENCE.out.fai,
  			ch0
  			)
  	versions = versions.mix(META_TO_BAMS.out.versions)


	bed = [[id:"nobed"],[]]
	if(params.bed!=null) {
		bed = [[id:"capture"],file(params.bed)];
		PREPARE_USER_BED(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.scatter_bed,
			bed
			)
		versions = versions.mix(PREPARE_USER_BED.out.versions)

		bed = PREPARE_USER_BED.out.bed
		}

	MANTA_SINGLE(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		META_TO_BAMS.out.bams,
		bed
		)
	}
runOnComplete(workflow)


