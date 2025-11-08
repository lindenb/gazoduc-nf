/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {BCFTOOLS_CALL as CALL } from '../../../modules/bcftools/call'
include {BCFTOOLS_MERGE        } from '../../../modules/bcftools/merge3'
include {BCFTOOLS_CONCAT       } from '../../../modules/bcftools/concat3'
include {isBlank               } from '../../../modules/utils/functions'


workflow BCFTOOLS_CALL {
	take:
		metadata
		fasta
		fai
		dict
		pedigree
		beds
		bams //[ meta, bam,bai]
	main:
		versions = Channel.empty()
		ch1 = bams
			.map{meta,bam,bai->
				if(!isBlank(meta.batch)) return [meta,bam,bai];
				return [meta.plus(batch:meta.id),bam,bai];
				}
			.map{
				if(isBlank(it[0].batch)) throw new IllegalArgumentException("Empty batch ??? ${it}");
				return it;
				}
			.map{meta,bam,bai->[[id:meta.batch],[bam,bai]]}
			.groupTuple()
			.map{meta,files->[ meta, files.flatten().sort() ]}
		
		ch2 = ch1.combine(beds.map{_meta,bed->bed}).view()

		CALL(
			fasta,
			fai,
			pedigree,
			[[id:"noploidy"],[]],
			ch2
			)
		versions = versions.mix(CALL.out.versions)

		/* group data by BATCH */
		group_by_batch = CALL.out.vcf
				.map{meta,vcf,idx,_bed->[meta,[vcf,idx]]}
				.groupTuple()
				.map{meta,files ->[meta, files.flatten().sort() ]}
				

		/** concat per batch, so we can have one vcf per batch/sample  */
		BCFTOOLS_CONCAT(group_by_batch)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

		// merge those having more than one sample
		BCFTOOLS_MERGE(
			BCFTOOLS_CONCAT.out.vcf
				.map{meta,vcf,tbi-> [ [id:"bcftools.call"],[vcf,tbi]]}
				.groupTuple()
				.map{meta,files ->[meta, files.flatten().sort() ]}
			)
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)
		vcf_out = BCFTOOLS_MERGE.out.vcf
		
	emit:
		versions
		vcf = vcf_out
	}


