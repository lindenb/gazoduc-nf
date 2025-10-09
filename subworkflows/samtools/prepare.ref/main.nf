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
include {SAMTOOLS_DICT         } from '../../../modules/samtools/dict'
include {SAMTOOLS_FAIDX        } from '../../../modules/samtools/faidx'
include {FAI2BED               } from '../../../modules/samtools/fai2bed'
include {SCATTER_TO_BED        } from '../../../subworkflows/gatk/scatterintervals2bed'
include {BEDTOOLS_COMPLEMENT   } from '../../../modules/bedtools/complement'


workflow PREPARE_REFERENCE {
take:
	meta
	fasta
main:
	versions=Channel.empty()

	SAMTOOLS_FAIDX(fasta)
	versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
	
	SAMTOOLS_DICT(fasta)
	versions = versions.mix(SAMTOOLS_DICT.out.versions)
	
	FAI2BED(SAMTOOLS_FAIDX.out.fai)
	versions = versions.mix(FAI2BED.out.versions)
	
	
	complement_bed =  Channel.of([[id:"nocomplement"],[]])
	if(meta.skip_scatter==null || meta.skip_scatter==false) {
		SCATTER_TO_BED(
			meta,
			fasta,
			SAMTOOLS_FAIDX.out.fai,
			SAMTOOLS_DICT.out.dict
			)
		versions = versions.mix(SCATTER_TO_BED.out.versions)
		scatter_bed = SCATTER_TO_BED.out.bed

		if(meta.skip_complement=null || meta.skip_complement==false) {
			BEDTOOLS_COMPLEMENT(
				SAMTOOLS_FAIDX.out.fai,
				SCATTER_TO_BED.out.bed
				)
			complement_bed = BEDTOOLS_COMPLEMENT.out.bed
			versions = versions.mix(BEDTOOLS_COMPLEMENT.out.versions)
			}
		}
	else
		{
		scatter_bed = Channel.of([[id:"noscatterbed"],[]])
		}
	

	

emit:
	versions
	fai = SAMTOOLS_FAIDX.out.fai
	dict = SAMTOOLS_DICT.out.dict
	bed = FAI2BED.out.bed
	scatter_bed
	complement_bed
}
