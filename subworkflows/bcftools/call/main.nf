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
include {BCFTOOLS_MERGE        } from '../../../modules/bcftools/merge'
include {BCFTOOLS_CONCAT       } from '../../../modules/bcftools/concat'



workflow BCFTOOLS_CALL {
	take:
		meta
		fasta
		fai
		dict
		pedigree
		beds
		bams //[ meta, bam,bai]
	main:
		versions = Channel.empty()
		ch1 = bams
			.map{
				if(it[0].containsKey("batch")) return it;
				return [it[0].plus(batch:it[0].sample),it[1],it[2]];
				}
			.map{[[id:it[0].batch],[it[1],it[2]]]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]}
		
		ch2 = ch1.combine(beds.map{it[1]})

		CALL(
			fasta,
			fai,
			pedigree,
			[[id:"noploidy"],[]],
			ch2
			)
		versions = versions.mix(CALL.out.versions)



		BCFTOOLS_MERGE(
			CALL.out.vcf
				.map{[[id:it[3].toRealPath().toString()/* bed */],[it[1],it[2]]]}
				.groupTuple()
				.map{[it[0],it[1].flatten(), [] /* no bed */]}
			)
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)

		BCFTOOLS_CONCAT(
			BCFTOOLS_MERGE.out.vcf
				.map{[[id:"call"],[it[1],it[2]]]}
				.groupTuple()
				.map{[it[0],it[1].flatten()]},
			[[id:"nobed"],[]]
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

	emit:
		versions
		vcf = BCFTOOLS_CONCAT.out.vcf
	}

