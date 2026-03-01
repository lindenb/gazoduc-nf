/*

Copyright (c) 2026 Pierre Lindenbaum

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
include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { DICT_TO_BED                              } from '../../modules/jvarkit/dict2bed'
include { runOnComplete                            } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { verify                                   } from '../../modules/utils/functions.nf'
include { VCF_TO_CONTIGS                           } from '../../subworkflows/bcftools/vcf2contigs'
include { DOWNLOAD_CHAIN                           } from '../../modules/ucsc/download.chain'
include {LIFTOVER_VCF                              } from '../../modules/gatk/liftovervcf'
include { MERGE_VCFS  as MERGE_VCFS_LIFTED         } from '../../modules/gatk/mergevcfs'
include { MERGE_VCFS  as MERGE_VCFS_FAIL           } from '../../modules/gatk/mergevcfs'

String toUcsc(String build) {
	if(build.equalsIgnoreCase("grch38")) return "hg38";
    if(build.equalsIgnoreCase("grch37")) return "hg19";
	return build;
	}

workflow {

		if(!workflow.stubRun) {
		validateParameters()
		}

	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	def metadata = [id:"liftover"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.vcf==null) {
			throw new IllegalArgumentException("undefined --vcf");
			}


  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/**
	* prepare VCF inputs :VCF may be composed of multiple VCFs (one per contig)
	*/
	VCF_INPUT(metadata.plus([
		path: params.vcf,
		arg_name: "vcf",
		require_index : true,
		required: true,
		unique : false
		]))
	versions = versions.mix(VCF_INPUT.out.versions)

	vcfs= VCF_INPUT.out.vcf
	DICT_TO_BED(vcfs.map{meta,vcf,_tbi->[meta,vcf]})
	versions = versions.mix(DICT_TO_BED.out.versions)
	
	ch1 = DICT_TO_BED.out.bed
		.splitCsv(header:true,sep:'\t')
		.map{meta,row->[meta,row.buildName]}
		.filter{_meta,build->!(isBlank(build) || build==".")}
		.map{meta,build->[meta,toUcsc(build)]}
		.unique()
		.groupTuple()
		.map{meta,builds->{
			verify(builds.size()==1,"multiple builds for ${meta}")
			return [meta,builds[0]];
			}}
		.map{meta,build->[meta,build]}
		.join(vcfs)
		.combine(PREPARE_ONE_REFERENCE.out.dict)
		.multiMap{meta1,build,vcf,_tbi,meta2,dict->
			source: [meta1.plus(ucsc_name:build),vcf]
			dest: [meta2,dict]
			}
	DOWNLOAD_CHAIN(ch1.source,ch1.dest)
	versions = versions.mix(DOWNLOAD_CHAIN.out.versions)
	

	VCF_TO_CONTIGS(metadata,vcfs)
	versions = versions.mix(VCF_TO_CONTIGS.out.versions)


	ch2 = DOWNLOAD_CHAIN.out.chain
		.map{meta1,_meta2,chain->[meta1,chain]}
		.combine(VCF_TO_CONTIGS.out.vcf)
		.filter{meta1,_chain,meta2,_vcf,_tbi->meta1.id==meta2.id}
		.multiMap{meta1,liftchain,meta2,vcffile,tbi->
			chain: [meta1,liftchain]
			vcf: [meta2,vcffile,tbi]
			}
	
	LIFTOVER_VCF(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		ch2.chain,
		ch2.vcf
		)
	versions = versions.mix(VCF_TO_CONTIGS.out.versions)

	MERGE_VCFS_LIFTED(
		LIFTOVER_VCF.out.vcf
			.map{meta,vcf,tbi->[[id:meta.id],[vcf,tbi]]}
			.groupTuple()
			.map{meta,files->[meta,files.flatten().sort()]}
		)
	versions = versions.mix(MERGE_VCFS_LIFTED.out.versions)

	MERGE_VCFS_FAIL(
		LIFTOVER_VCF.out.fail
			.map{meta,vcf,tbi->[[id:meta.id],[vcf,tbi]]}
			.groupTuple()
			.map{meta,files->[meta,files.flatten().sort()]}
		)
	versions = versions.mix(MERGE_VCFS_FAIL.out.versions)
	}

runOnComplete(workflow)
