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
include { WGSELECT                                 } from '../../modules/wgselect'
include { SNPEFF_DOWNLOAD                          } from '../../modules/snpeff/download'
include { BCFTOOLS_CONCAT                          } from '../../modules/bcftools/concat3'
include { MULTIQC as MQC1                          } from '../../modules/multiqc'
include { MULTIQC as MQC2                          } from '../../modules/multiqc'
include { JVARKIT_MULTIQCPOSTPROC                  } from '../../modules/jvarkit/multiqcpostproc'

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
	if(params.fasta==null) {
		log.warn("undefined --fasta")
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
		log.warn("undefined --cases")
		}
	else{
		cases_ch = [[id:"cases"],file(params.cases)]
		}

	if(params.controls==null) {
		controls_ch = [[id:"controls"],[]]
		log.warn("undefined --controls")
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
	bed = bed.first()
	
	black_list_bed_ch = [[id:"noblacklist"],[]]
	if(parseBoolean(params.with_encode_blacklist)) {
		/***************************************************
		*
		* Download encode black list
		*
		*/
		ENCODE_BLACKLIST(PREPARE_ONE_REFERENCE.out.dict)
		versions = versions.mix(ENCODE_BLACKLIST.out.versions)
		black_list_bed_ch = ENCODE_BLACKLIST.out.bed
		
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
	
	vcf = VCF_INPUT.out.vcf
		.map{meta,vcf,tbi->[meta.plus(id:vcf.toRealPath().toString().md5()),vcf,tbi]}
		

	VCF_TO_BED(vcf)
	versions = versions.mix(VCF_INPUT.out.versions)

	vcf = VCF_TO_BED.out.bed
		.splitCsv(header:false,sep:'\t')
		.map{meta,row->[meta,row[0]]}
		.combine(vcf)
		.filter{meta1,_contig,meta2,_vcf,_tbi->meta1.id==meta2.id}
		.map{meta1,ctg,_meta2,vcffile,tbi->[meta1.plus(contig:ctg),vcffile,tbi]}
		
	
	SPLIT_N_VARIANTS( bed, vcf )
	versions = versions.mix(VCF_INPUT.out.versions)

	vcfs = SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
            .combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
            .filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
            .map{meta1,vcf,meta2,tbi->[meta1.plus(id:vcf.name.md5()),vcf,tbi]}
			
	
	if(params.gnomad==null) {
		gnomad_vcf = [[id:"nognomad"],[],[]]
		}
	else
		{
		gnomad_vcf = [[id:"gnomad"],file("${params.gnomad}"),file("${params.gnomad}.tbi")]
		}
	
	SNPEFF_DOWNLOAD(
		PREPARE_ONE_REFERENCE.out.fai
			.map{meta,_fai->[meta,file(params.snpeff_database_directory)]}
		)
	versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)


	if(params.gatk_hard_filters_arguments==null) {
		apply_hard_filters_arguments = [[id:"nogatkhardfilters"],[]]
		}
	else
		{
		apply_hard_filters_arguments = [[id:"gatkhardfilters"],file(params.gatk_hard_filters_arguments)]
		}

	cadd_ch = [[id:"nocadd"],[],[]]
	
	if(params.mappability_bigwig==null) {
		mappability = [[id:"nomappability"],[]]
		}
	else
		{
		mappability = [[id:"mappability"],file("${params.mappability_bigwig}")]
		}

	WGSELECT(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		gnomad_vcf,
		SNPEFF_DOWNLOAD.out.database,
		cadd_ch,
		mappability,
		black_list_bed_ch,
		apply_hard_filters_arguments,
		cases_ch,
		controls_ch,
		vcfs.map{meta,vcffile,tbi->[meta,vcffile,tbi,[]]}
		)
	versions = versions.mix(WGSELECT.out.versions)

	BCFTOOLS_CONCAT(
		WGSELECT.out.vcf
			.map{_meta,vcf,tbi->[vcf,tbi]}
			.flatMap()
			.collect()
			.map{files->[[id:"wgselect"],files.sort()]}
		)
	versions = versions.mix(BCFTOOLS_CONCAT.out.versions)


	def gtf    = [[id:"gtf"], file(params.gtf), file(params.gtf+".tbi") ]
	def gff3    = [[id:"gff3"], file(params.gff3), file(params.gff3+".tbi") ]

	MAKE_PEDIGREE(cases_ch,controls_ch)
	versions = versions.mix(MAKE_PEDIGREE.out.versions)

	VCF_STATS_AFTER(
		metadata.plus(
			with_bcftools_stats  : true,
			gatk_denovo : false,
			ad_ratio : true
		),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		Channel.of(gtf),
		Channel.of(gff3),
		MAKE_PEDIGREE.out.pedigree.ifEmpty([[id:"nopedigree"],[]]),
		MAKE_PEDIGREE.out.sample2pop.ifEmpty([[id:"nogroup2sample"],[]]),
		bed,
		BCFTOOLS_CONCAT.out.vcf
		)
	versions = versions.mix(VCF_STATS_AFTER.out.versions)
	multiqc = multiqc.mix(VCF_STATS_AFTER.out.multiqc)



	MQC1(
		[[id:"no_mqc_config"],[]],
		multiqc.map{_meta,f->f}.collect().map{files->[metadata,files.sort()]}
		)
	
	JVARKIT_MULTIQCPOSTPROC(
		MAKE_PEDIGREE.out.sample2pop,
		MQC1.out.datadir
		)
	
	MQC2(
		[[id:"no_mqc_config"],[]],
		JVARKIT_MULTIQCPOSTPROC.out.json.map{_meta,f->f}.collect().map{files->[metadata.plus(id:"perpop"),files.sort()]}
		)

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

process MAKE_PEDIGREE {
label "process_single"
input:
	tuple val(meta1),path(cases)
	tuple val(meta2),path(controls)
output:
	tuple val(meta1),path("*.ped"),optional:true,emit:pedigree
	tuple val(meta1),path("*.tsv"),optional:true,emit:sample2pop
	path("versions.yml"),emit:versions
script:
"""
mkdir -p TMP
if ${(cases?true:false) && (controls?true:false)}
then
	sort "${cases}" |uniq > TMP/cases.txt
	sort "${controls}" |uniq > TMP/controls.txt
	awk -F '\t' '{printf("%s\t%s\t0\t0\t0\tcase\\n",\$1,\$1);}' TMP/cases.txt >  TMP/pedigree.ped
	awk -F '\t' '{printf("%s\t%s\t0\t0\t0\tcontrol\\n",\$1,\$1);}' TMP/controls.txt >> TMP/pedigree.ped
	
	
	awk -F '\t' '{printf("%s\tcase\\n",\$1);}' TMP/cases.txt >  TMP/sample2pop.tsv
	awk -F '\t' '{printf("%s\tcontrol\\n",\$1);}' TMP/controls.txt >> TMP/sample2pop.tsv

	mv TMP/pedigree.ped ./
	mv TMP/sample2pop.tsv ./
fi
touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}