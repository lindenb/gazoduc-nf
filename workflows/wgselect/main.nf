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


include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { parseBoolean                             } from '../../modules/utils/functions.nf'
include { flatMapByIndex                           } from '../../modules/utils/functions.nf'
//include { WGSELECT_02                              } from '../../subworkflows/wgselect/wgselect2'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect' 
include { PREPARE_USER_BED                         } from '../../subworkflows/bedtools/prepare.user.bed'
include { ENCODE_BLACKLIST                         } from '../../modules/encode/blacklist' 
include { BEDTOOLS_SUBTRACT                        } from '../../modules/bedtools/subtract' 
include { SPLIT_N_VARIANTS                         } from '../../modules/jvarkit/splitnvariants'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { VCF_TO_BED                               } from '../../modules/bcftools/vcf2bed'
include { VCF_STATS as VCF_STATS_BEFORE            } from '../../subworkflows/vcfstats'
include { VCF_STATS as VCF_STATS_AFTER             } from '../../subworkflows/vcfstats'


workflow {
	
	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}

	
		
	versions = Channel.empty()
	multiqc = Channel.empty()
	metadata = [id:"wgselect"]
	
	if(params.vcf==null) {
		log.warn("undefined --vcf")
		exit -1
		}

	
	if(params.vcf==null) {
		log.warn("undefined --vcf")
		exit -1
		}


  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata,
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	if(params.cases==null) {
		cases_ch = [[id:"no_cases"],[]]
		}
	else{
		cases_ch = [[id:"cases"],file(params.cases)]
		}

	if(params.controls==null) {
		controls_ch = [[id:"controls"],[]]
		}
	else{
		controls_ch = [[id:"controls"],file(params.controls)]
		}



	/***************************************************
	*
	*  MAKE BED
	*
	*/
	if(params.bed==null) {
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
		}
	else
		{
		BEDTOOLS_INTERSECT(
			PREPARE_ONE_REFERENCE.out.fai,
			 Channel.of(params.bed)
			 	.map{bed->file(bed)}
				.map{bed->[[id:bed.baseName],bed]}
				.combine(PREPARE_ONE_REFERENCE.out.scatter_bed)
				.map{meta1,bed1,meta2,bed2->[meta1,bed1,bed2]}
			)
		versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)
		bed = BEDTOOLS_INTERSECT.out.bed
		}
	
	if(parseBoolean(params.with_encode_blacklist)) {
		/***************************************************
		*
		* Download encode black list
		*
		*/
		ENCODE_BLACKLIST(PREPARE_ONE_REFERENCE.out.dict)
		versions = versions.mix(ENCODE_BLACKLIST.out.versions)
		BEDTOOLS_SUBTRACT(bed.combine(ENCODE_BLACKLIST.out.bed).map{meta1,bed1,_meta2,bed2->[meta1,bed1,bed2]})
		versions = versions.mix(BEDTOOLS_SUBTRACT.out.versions)
		bed = BEDTOOLS_SUBTRACT.out.bed.first()
		}

	
	/***************************************************
	*
	* VCF input
	*
	*/
	VCF_INPUT(metadata.plus([
		path: params.vcf,
		arg_name: "vcf",
		require_index : true,
		required: true,
		unique : false
		]))
	versions = versions.mix(VCF_INPUT.out.versions)
	
	SPLIT_N_VARIANTS(
		bed,
		VCF_INPUT.out.vcf
		)
	versions = versions.mix(VCF_INPUT.out.versions)

	vcfs = SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
            .combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
            .filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
            .map{meta1,vcf,meta2,tbi->[meta1.plus(id:vcf.name.md5()),vcf,tbi]}



	

	//VCF_TO_BED.out.bed.
	
	/*
	WGSELECT_02(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		Channel.fromPath(params.vcf),
		file(params.samples),
		file(params.bed)
		)
	versions = versions.mix(WGSELECT_02.out.versions)
	multiqc = versions.mix(WGSELECT_02.out.multiqc)
	*/
	}

