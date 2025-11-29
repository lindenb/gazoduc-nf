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
include { parseBoolean                } from '../../modules/utils/functions.nf'
include { assertKeyExistsAndNotEmpty  } from '../../modules/utils/functions.nf'
include { META_TO_BAMS                } from '../../subworkflows/samtools/meta2bams1'
include { MANTA                       } from '../../subworkflows/manta' 
include { PREPARE_ONE_REFERENCE       } from '../../subworkflows/samtools/prepare.one.ref'
include { READ_SAMPLESHEET            } from '../../subworkflows/nf/read_samplesheet'
include { PREPARE_USER_BED            } from '../../subworkflows/bedtools/prepare.user.bed'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}
	multiqc = Channel.empty()
	versions = Channel.empty()
	def metadata = [id:"manta"]
		

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	
	READ_SAMPLESHEET(
		[arg_name:"samplesheet"],
		params.samplesheet
		)
	versions = versions.mix(READ_SAMPLESHEET.out.versions)
	
	META_TO_BAMS(
		metadata ,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		READ_SAMPLESHEET.out.samplesheet
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

		bed = PREPARE_USER_BED.out.bed.first()
		}

	MANTA(
		metadata.plus(with_truvari: parseBoolean(params.with_truvari)),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		META_TO_BAMS.out.bams,
		bed
		)
	versions = versions.mix(MANTA.out.versions)
	multiqc = multiqc.mix(MANTA.out.versions)
	}
runOnComplete(workflow)


