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
include {BCFTOOLS_CONCAT                      } from '../../../modules/bcftools/concat2'
include {BCFTOOLS_MERGE                       } from '../../../modules/bcftools/merge'

workflow HAPLOTYPECALLER_DIRECT {
	take:
		meta
		fasta
		fai
		dict
		dbsnp
		pedigree
		beds
		bams
	main:

		version_ch = Channel.empty()



		ch1 = bams.map{meta,bam,bai->["${meta.group?:meta.id}", [bam, bai]]}
			.groupTuple()
			.map{group,bam_files->[[id:group],bam_files.flatten().sort()]}
			.combine(beds.map{meta,bed->bed})
			
		
		HAPCALLER(
			fasta,
			fai,
			dict,
			dbsnp,
			pedigree,
			ch1
			)
		version_ch = version_ch.mix(HAPCALLER.out.versions)


		
		to_merge = HAPCALLER.out.vcf
				.map{meta,vcf,tbi,bed->[bed.toRealPath().toString(),[vcf,tbi]]}
				.groupTuple()
				.map{name,vcf_files->[[id:name.md5()],vcf_files.flatten().sort()]}
				
				
		BCFTOOLS_MERGE(
			[[id:"nobed"],[]],
			to_merge
			)
			
		version_ch = version_ch.mix(BCFTOOLS_MERGE.out.versions)

		BCFTOOLS_CONCAT(
			[[id:"nobed"],[]],
			BCFTOOLS_MERGE.out.vcf
				.map{meta,vcf,tbi->[vcf,tbi]}
				.flatMap()
				.collect()
				.map{[
					[id:meta.id],
					it
					]}
			)
		version_ch = version_ch.mix(BCFTOOLS_CONCAT.out.versions)
		
		vcf_out = BCFTOOLS_CONCAT.out.vcf.map{meta,vcf,tbi,bed->[meta,vcf,tbi,[]]}
		TODO AFTER_CONCAT
	emit:
		versions = version_ch
		vcf = vcf_out
	}

