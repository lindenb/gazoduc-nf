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

include { PREPARE_ONE_REFERENCE                } from '../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                            } from '../../subworkflows/nf/vcf_input'
include { runOnComplete                        } from '../../modules/utils/functions.nf'
include { parseBoolean                         } from '../../modules/utils/functions.nf'
include { isBlank                              } from '../../modules/utils/functions.nf'
include { TRUVARI                              } from '../../subworkflows/truvari'
include { BCFTOOLS_INDEX as  BCFTOOLS_INDEX1   } from '../../modules/bcftools/index'
include { BCFTOOLS_INDEX as  BCFTOOLS_INDEX2   } from '../../modules/bcftools/index'
include { BCFTOOLS_CONCAT                      } from '../../modules/bcftools/concat3'
include { JVARKIT_VCF_SET_DICTIONARY           } from '../../modules/jvarkit/vcfsetdict'
include { MULTIQC                              } from '../../modules/multiqc'
include { COMPILE_VERSIONS                     } from '../../modules/versions'
include { JVARKIT_VCFSTATS                     } from '../../modules/jvarkit/vcfstats'
include { VCF_TO_CONTIGS                       } from '../../subworkflows/bcftools/vcf2contigs'
include { GTF_INPUT                            } from '../../subworkflows/nf/gtf_input'
include { GTF_ANNOTATION                       } from '../../modules/gtf/annot1'
include { DOWNLOAD_GNOMAD_SV                   } from '../../modules/gnomad_sv/download.vcf'
include { JVARKIT_VCFGNOMADSV                  } from '../../modules/jvarkit/vcfgnomadsv'
include { BCFTOOLS_VIEW                        } from '../../modules/bcftools/view'
include { BCFTOOLS_QUERY_SAMPLES               } from '../../modules/bcftools/query.samples'
workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()


	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}

	def metadata = [id:"truvari"]
		

	PREPARE_ONE_REFERENCE(
		metadata.plus(skip_scatter:true),
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	
	VCF_INPUT(
		metadata.plus(
			path: "${params.samplesheet}",
			arg_name : "samplesheet",
			require_index : false,
			required: true,
			unique : false
			)
		)
	versions = versions.mix(VCF_INPUT.out.versions)
	vcf_ch = VCF_INPUT.out.vcf

	/*
	* check there is no duplicate sample in the vcf files
	*/
	if(parseBoolean(params.with_check_no_dup_sample)) {
	/** check uniq samples */
	BCFTOOLS_QUERY_SAMPLES(vcf_ch)
	versions = versions.mix(BCFTOOLS_QUERY_SAMPLES.out.versions)
	BCFTOOLS_QUERY_SAMPLES.out.samples.splitText()
		.map{meta,sn->[sn.trim(),meta.id]}
		.filter{sn,_meta->!isBlank(sn)}
		.groupTuple()
		.filter{_sn,array->array.size()>1}
		.subscribe{sn,array->
			def msg = "duplicate sample ${sn} in ${array}";
			log.error(msg);
			throw new IllegalArgumentException(msg);
			}
	}
	/**
	 *
	 * Mixing different versions of the same build ?
	 */
	if(parseBoolean(params.with_vcf_set_dict)) {
		JVARKIT_VCF_SET_DICTIONARY(
			PREPARE_ONE_REFERENCE.out.dict,
			vcf_ch.map{meta,f->[meta,f,[]]}
			)
		versions = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions)
		vcf_ch = JVARKIT_VCF_SET_DICTIONARY.out.vcf.map{meta,f,tbi->[meta,f]}
		}

	if(params.bed!=null) {
		BCFTOOLS_VIEW(
			[[id:"capture"],file(params.bed)],
			[[id:"nosamples"],[]],
			vcf_ch.map{meta,vcf->[meta,vcf,[]]}
			)
		versions = versions.mix(BCFTOOLS_VIEW.out.versions)
		vcf_ch = BCFTOOLS_VIEW.out.vcf
		}

	BCFTOOLS_INDEX1(vcf_ch)
	versions = versions.mix(BCFTOOLS_INDEX1.out.versions)
	vcf_ch = BCFTOOLS_INDEX1.out.vcf 

	TRUVARI(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		vcf_ch
		)
	versions = versions.mix(TRUVARI.out.versions)
	multiqc = multiqc.mix(TRUVARI.out.multiqc)
	vcf_ch = TRUVARI.out.vcf

	if(parseBoolean(params.with_annotation)) {

		GTF_INPUT(
			metadata.plus(
				arg_name:"gtf",
				require_index:false,
				download:true,
				path: params.gtf
				),
			PREPARE_ONE_REFERENCE.out.dict
			)
		versions = versions.mix(GTF_INPUT.out.versions)

		GTF_ANNOTATION(
			PREPARE_ONE_REFERENCE.out.fai,
			GTF_INPUT.out.gtf,
			[[id:"nogeneofinterest"],[]],
			vcf_ch.map{meta,vcf,tbi->[meta,vcf]}
			)
		versions = versions.mix(GTF_ANNOTATION.out.versions)
		vcf_ch = GTF_ANNOTATION.out.vcf
		
		DOWNLOAD_GNOMAD_SV( PREPARE_ONE_REFERENCE.out.dict )
		versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)
		
		JVARKIT_VCFGNOMADSV(
				DOWNLOAD_GNOMAD_SV.out.vcf,
				vcf_ch
				)
		versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)
		vcf_ch = JVARKIT_VCFGNOMADSV.out.vcf

		BCFTOOLS_INDEX2(vcf_ch)
		versions = versions.mix(BCFTOOLS_INDEX2.out.versions)
		vcf_ch = BCFTOOLS_INDEX2.out.vcf
		}


	/*run statistics over the vcf */
	JVARKIT_VCFSTATS(
		[[id:"nosample2pop"],[]],
		vcf_ch
		)
	versions = versions.mix(JVARKIT_VCFSTATS.out.versions)
    multiqc = multiqc.mix(JVARKIT_VCFSTATS.out.json)

	if(parseBoolean(params.with_multiqc)) {
		/**
		* At the end, make QC, make versions, zip results
		*/
		COMPILE_VERSIONS(versions.collect().map{it.sort()})

		MULTIQC(
			[[id:"noconfig"],[]],
			multiqc
				.map{_meta,f->f}
				.mix(COMPILE_VERSIONS.out.multiqc)
				.collect()
				.map{files->[metadata,files.sort()]}
			)
		}
	}

runOnComplete(workflow)
