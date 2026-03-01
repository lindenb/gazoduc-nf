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
include { GATK_GENOTYPE_VCF                       } from '../../modules/gatk/genotypevcf'
include { JVARKIT_VCFFILTERJDK                    } from '../../modules/jvarkit/vcffilterjdk'
include { BCFTOOLS_MERGE                          } from '../../modules/bcftools/merge3'
include { BCFTOOLS_INDEX                          } from '../../modules/bcftools/index'

workflow VS_CONTROL {
take:
    metadata
    fasta
    fai
    dict
    jvarkit_filter // [meta,filter]
    user_vcf // [meta,vcf,tbi]
    bam // [meta,bam,bai]
main:
    versions = Channel.empty()
    multiqc =  Channel.empty()

    GATK_GENOTYPE_VCF(
        fasta,
        fai,
        dict,
        [[id:"nodbsnp"],[],[]],
        user_vcf,
        bam.map{meta,bam,bai->[meta,[bam,bai]]}
        )
    versions = versions.mix(GATK_GENOTYPE_VCF.out.versions)

    BCFTOOLS_MERGE(GATK_GENOTYPE_VCF.out.vcf
        .combine(user_vcf)
        .map{meta1,vcf1,tbi1,_meta2,vcf2,tbi2->[meta1,[vcf1,tbi1,vcf2,tbi2]]})
    versions = versions.mix(BCFTOOLS_MERGE.out.versions)

	JVARKIT_VCFFILTERJDK(
		jvarkit_filter,
		[[id:"no_ped"],[]],
		BCFTOOLS_MERGE.out.vcf
		)
	versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)
	
    BCFTOOLS_INDEX(JVARKIT_VCFFILTERJDK.out.vcf)
    versions = versions.mix(BCFTOOLS_INDEX.out.versions)

emit:
    vcf = BCFTOOLS_INDEX.out.vcf
    versions
    multiqc
}