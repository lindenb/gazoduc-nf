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

include { validateParameters         } from 'plugin/nf-schema'
include { paramsHelp                 } from 'plugin/nf-schema'
include { paramsSummaryLog           } from 'plugin/nf-schema'
include { samplesheetToList          } from 'plugin/nf-schema'
include { INDEXCOV                   } from '../../subworkflows/indexcov/simple'
include { JVARKIT_BAM_RENAME_CONTIGS } from '../../modules/jvarkit/bamrenamechr'
include { runOnComplete              } from '../../modules/utils/functions.nf'
include { COMPILE_VERSIONS           } from '../../modules/versions'
include { READ_SAMPLESHEET           } from '../../subworkflows/nf/read_samplesheet'
include { META_TO_BAMS               } from '../../subworkflows/samtools/meta2bams2'
include { META_TO_PED                } from '../../subworkflows/pedigree/meta2ped'
include {PREPARE_ONE_REFERENCE       } from '../../subworkflows/samtools/prepare.one.ref'




workflow {


	validateParameters()


	if( params.help ) {
		log.info(paramsHelp())
		exit 0
	}  else {
	// Print summary of supplied parameters
	log.info paramsSummaryLog(workflow)
	}


	/* no fastq samplesheet */
	if(params.samplesheet==null) {
		log.error("--samplesheet undefined")
		exit -1
		}
	/* no fastq samplesheet */
	if(params.fasta==null) {
		log.error("--fasta undefined")
		exit -1
		}




	versions = Channel.empty()
 	def metadata = [ id: "indexcov"   ]
	

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)



	samplesheet0_ch = READ_SAMPLESHEET(
		[arg_name:"samplesheet"],
		params.samplesheet
		).samplesheet
	versions = versions.mix(READ_SAMPLESHEET.out.versions)


	if(params.pedigree!=null) {
		pedigree = Channel.of([[id:"ped"],file(params.pedigree)])
		}
	else
		{
		META_TO_PED(metadata, samplesheet0_ch.map{it[0]})
		versions = versions.mix(META_TO_PED.out.versions)
		pedigree = META_TO_PED.out.pedigree_gatk
		}
	
	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)


  	bams = META_TO_BAMS.out.bams
			.branch {meta,bam,bai,fasta,fai,dict->
				ok_ref: fasta.toRealPath().toString().equals(params.fasta)
				bad_ref: true
				}
			//[[],]}
   
	
	
	JVARKIT_BAM_RENAME_CONTIGS(
		PREPARE_ONE_REFERENCE.out.dict,
		[[id:"nobed"],[]],
		bams.bad_ref
		)
	versions =versions.mix(JVARKIT_BAM_RENAME_CONTIGS.out.versions)

 	

	INDEXCOV(
		metadata.plus([
			batch_size:  (params.batch_size as int),
			radius: (params.radius as int)
			]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		PREPARE_ONE_REFERENCE.out.complement_bed,
		pedigree,
		bams.ok_ref.map{meta,bam,bai,_fasta,_fai,_dict->[[id:meta.id],bam,bai]}
			.mix(JVARKIT_BAM_RENAME_CONTIGS.out.bam.map{[it[0].plus([force_rebuild:false]),it[1],it[2]]})
		)
	versions = versions.mix(INDEXCOV.out.versions)

	COMPILE_VERSIONS(versions.collect())
	}


runOnComplete(workflow)
