
/*

THIS FILE WAS GENERATED DO NOT EDIT !!!

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



include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { ENCODE_BLACKLIST                         } from '../../modules/encode/blacklist' 
include { BEDTOOLS_SUBTRACT                        } from '../../modules/bedtools/subtract' 
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect' 
include { BED_CLUSTER                              } from '../../modules/jvarkit/bedcluster'
include { HAPLOTYPECALLER                          } from '../../subworkflows/gatk/haplotypecaller'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include { JVARKIT_VCFFILTERJDK                     } from '../../modules/jvarkit/vcffilterjdk'
include { BCFTOOLS_NORM                            } from '../../modules/bcftools/norm'
include { BCFTOOLS_INDEX                           } from '../../modules/bcftools/index'
include { MERGE_VCFS                               } from '../../modules/gatk/mergevcfs'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { GATK_GENOTYPE_VCF_RECURSIVE              } from '../../modules/gatk/genotypevcf.recursive'


workflow {
  versions = Channel.empty()
  multiqc = Channel.empty()


	def metadata = [id:"ultrares"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	
    if(params.jvarkit_filter==null) {
			throw new IllegalArgumentException("undefined --jvarkit_filter");
			}
	 if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
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



  /***************************************************
   *
   *  READ BAM SAMPLESSHEET
   *
   */
	READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        )
	versions = versions.mix(READ_SAMPLESHEET.out.versions)

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)

	if(params.vcf==null) {

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


		BEDTOOLS_MAKEWINDOWS( bed )
		versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)


		/* if it's an exome , group the small genome together in BED */
		BED_CLUSTER(
			PREPARE_ONE_REFERENCE.out.dict,
			BEDTOOLS_MAKEWINDOWS.out.bed
			)
		versions = versions.mix(BED_CLUSTER.out.versions)
		bed = BED_CLUSTER.out.bed
			.map{_meta,bed->bed}
			.map{it instanceof List?it:[it]}
			.flatMap()
			.map{bed->[[id:bed.baseName],bed]}
		

		/***************************************************
		*
		* Run Hapcaller on CASES
		*
		*/
	HAPLOTYPECALLER(
		metadata.plus(
			gvcf_merge_method : "combinegvcfs",
			with_split_bed: false
			),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		[[id:"oneref"],[]],
		[[id:"nodbsnp"],[],[]],
		[[id:"noped"],[]],
		bed,
		META_TO_BAMS.out.bams.filter{meta,_bam,_bai->meta.status=="case"}
		)
	versions = versions.mix(HAPLOTYPECALLER.out.versions)
	multiqc = multiqc.mix(HAPLOTYPECALLER.out.multiqc)
	vcf = 	HAPLOTYPECALLER.out.vcf_bed_ch.map{meta,vcf,tbi,_bed->[meta,vcf,tbi]}
	}
else
	{
	vcf = Channel.of([[id:"user_vcf"],file(params.vcf),file(params.vcf+".tbi")])
	}
  
 
   BCFTOOLS_NORM(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		vcf
		)
	versions = versions.mix(BCFTOOLS_NORM.out.versions)
	vcf = BCFTOOLS_NORM.out.vcf

   jvarkit_filter = [[id:"jvarkit_filter"],file(params.jvarkit_filter)]
	JVARKIT_VCFFILTERJDK(
		jvarkit_filter,
		[[id:"no_ped"],[]],
		vcf
		)
	versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)
	
    BCFTOOLS_INDEX(JVARKIT_VCFFILTERJDK.out.vcf)
    versions = versions.mix(BCFTOOLS_INDEX.out.versions)

	vcf = BCFTOOLS_INDEX.out.vcf

	GATK_GENOTYPE_VCF_RECURSIVE(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		[[id:"no_dbsnp"],[],[]],
		jvarkit_filter,
		vcf,
		META_TO_BAMS.out.bams
			.filter{meta,_bam,_bai->meta.status!="case"}
			.map{_meta,bam,bai->[bam,bai]}
			.flatMap()
			.collect()
			.map{files->[metadata,files.sort()]}
		)
	versions = versions.mix(GATK_GENOTYPE_VCF_RECURSIVE.out.versions)

  MERGE_VCFS(
	GATK_GENOTYPE_VCF_RECURSIVE.out.vcf.map{meta,vcf,tbi->[vcf,tbi]}
		.flatMap()
		.collect()
		.map{files->[metadata,files.sort()]}
  	)
  versions = versions.mix(MERGE_VCFS.out.versions)
}
