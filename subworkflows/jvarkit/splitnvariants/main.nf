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
include {SPLIT_N_VARIANTS as SPLIT2VARIANTS } from '../../../modules/jvarkit/splitnvariants/main.nf'
include {VCF_TO_BED                         }  from '../../../modules/bcftools/vcf2bed/main.nf'
include {BED_CLUSTER                        }  from '../../../modules/jvarkit/bedcluster/main.nf'
include {BEDTOOLS_INTERSECT                 }  from '../../../modules/bedtools/intersect/main.nf'


/**
 * call jvarkit vcfsplitnvariants one or more vcf
 *
 */
workflow SPLIT_N_VARIANTS {
	take:
		meta
		fasta
		fai
		dict
		vcf
		optional_bed
	main:
		versions = Channel.empty()
		bed_ch = VCF_TO_BED(vcf)
		versions.mix(VCF_TO_BED.out.versions)

		bed_ch = VCF_TO_BED.out.bed

		if(optional_bed[1]) {
			BEDTOOLS_INTERSECT(fai, bed_ch.map{[it[0],it[1],optional_bed[1]]})
			versions.mix(BEDTOOLS_INTERSECT.out.versions)
			bed_ch = BEDTOOLS_INTERSECT.out.bed
		}

		BED_CLUSTER(fasta, fai, dict, VCF_TO_BED.out.bed)
		versions.mix(BED_CLUSTER.out.versions)

		ch1 = BED_CLUSTER.out.bed
			.join(vcf)
			.view()
		
		SPLIT2VARIANTS(ch1)

		vcf = SPLIT2VARIANTS.out.vcf
			.join(SPLIT2VARIANTS.out.tbi)
			.view()
/*
		bed_ch = VCF_TO_BED(vcf,optional_bed)
		versions.mix(VCF_TO_BED.out.versions)


		vcf = VCF_TO_BED.out.vcf.join(VCF_TO_BED.out.vcf.tbi)
		versions.mix(VCF_TO_BED.out.versions)*/
	emit:
		vcf
		versions = versions
	}


