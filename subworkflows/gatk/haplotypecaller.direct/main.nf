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
//include {makeKey                } from '../../../modules/utils/functions.nf'
include {HAPLOTYPECALLER_DIRECT as HAPCALLER  } from '../../../modules/gatk/hapcallerdirect'
include {BCFTOOLS_CONCAT                      } from '../../../modules/bcftools/concat3'
include {BCFTOOLS_MERGE                       } from '../../../modules/bcftools/merge3'
include {isBlank                              } from '../../../modules/utils/functions'
include {parseBoolean                         } from '../../../modules/utils/functions'

workflow HAPLOTYPECALLER_DIRECT {
	take:
		metadata
		fasta
		fai
		dict
		dbsnp
		pedigree
		beds
		bams
	main:

		version_ch = Channel.empty()
		/* checl all beds have an unique ID */
		beds.map{meta,bed->[meta.id,bed]}
			.groupTuple()
			.filter{it[1].size()!=1}
			.map{throw new IllegalArgumentException("In HAPLOTYPECALLER_DIRECT bed should have a unique id"); return it;}


		bams.filter{meta,_bam,_bai->isBlank(meta.batch)}
			.take(1)
			.map{log.warn("in HAPLOTYPECALLER_DIRECT : meta is missing a 'batch' property");return 0;}


		/* group all the bams, make a (meta,bams,bed] */
		ch1 = bams.map{meta,bam,bai->["${meta.batch?:meta.id}", [bam, bai]]}
			.groupTuple()
			.map{group_id,bam_files->[[batch:group_id],bam_files.flatten().sort()]}
			.combine(beds)
			.map{meta1,bam_files,meta2,bed->[meta1.plus(id:meta2.id),bam_files,bed]} //give them the bed.id
		
		HAPCALLER(
			fasta,
			fai,
			dict,
			dbsnp,
			pedigree,
			ch1
			)
		version_ch = version_ch.mix(HAPCALLER.out.versions)


		to_concat = HAPCALLER.out.vcf
				.map{meta,vcf,tbi,bed->[meta.batch,[vcf,tbi]]} // group by bedid
				.groupTuple()
				.map{batch_id,vcf_files->[[id:batch_id],vcf_files.flatten().sort()]}

		// concat per GROUP
		BCFTOOLS_CONCAT(to_concat)
		version_ch = version_ch.mix(BCFTOOLS_CONCAT.out.versions)
		

		if(metadata.with_merge==null || parseBoolean(metadata.with_merge)) {
			BCFTOOLS_MERGE(
				BCFTOOLS_CONCAT.out.vcf
					.map{_meta,vcf,idx->[vcf,idx]}
					.collect()
					.map{vcf_files->[[id:(isBlank(metadata.id)?"hapcaller.direct":metadata.id)],vcf_files.flatten().sort()]}
				)
				
			version_ch = version_ch.mix(BCFTOOLS_MERGE.out.versions)
			vcf_out = BCFTOOLS_MERGE.out.vcf
			}
		else
			{
			vcf_out = Channel.empty()
			}


	emit:
		versions = version_ch
		per_group_vcf = BCFTOOLS_CONCAT.out.vcf
		vcf = BCFTOOLS_MERGE.out.vcf
	}

