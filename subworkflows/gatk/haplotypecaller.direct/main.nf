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
include {isBlank                              } from '../../../modules/utils/functions'
include {parseBoolean                         } from '../../../modules/utils/functions'
include {SIMPLE_CALLING                       } from '../../../subworkflows/simple.calling'

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
		
		SIMPLE_CALLING(
			metadata.plus(method:"gatk"),
			fasta,
			fai,
			dict,
			dbsnp,
			pedigree,
			beds,
			bams
			)
		versions = versions.mix(SIMPLE_CALLING.out.versions)

	emit:
		versions = version_ch
		per_group_vcf = SIMPLE_CALLING.out.per_group_vcf
		vcf = SIMPLE_CALLING.out.vcf
	}

