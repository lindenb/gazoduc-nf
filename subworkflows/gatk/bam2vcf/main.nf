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

include {GATK_BAM2VCF as BAM2VCF                    } from '../../../modules/gatk/bam2vcf'
include {BCFTOOLS_CONCAT                            } from '../../../modules/bcftools/concat3'


workflow GATK_BAM2VCF {
 take:
    metadata
    fasta
    fai
    dict
    dbsnp
    pedigree
    references //[meta, [ref files fa fai dict...]] all known reference
    beds // [meta,bed]
    bams // [meta,bam,bai]
main:
    versions = Channel.empty()

    /* check all beds have an unique ID */
	beds.map{meta,bed->[meta.id,bed]}
		.groupTuple()
		.filter{it[1].size()!=1}
		.map{throw new IllegalArgumentException("In GATK_BAM2VCF bed should have a unique id"); return it;}

    BAM2VCF(
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        bams.combine(beds)
            .map{meta1,bam,bai,meta2,bed->[
                meta2, /* bed.meta */
                [bam,bai],/*bam and bai */
                bed
                ]}
            .groupTuple()
            .map{meta,bam_files,beds->[meta,bam_files.flatten().sort(), beds.sort()[0]]}
        )
    versions = versions.mix(BAM2VCF.out.versions)


    BCFTOOLS_CONCAT(
            BAM2VCF.out.vcf
            .map{meta,vcf,idx,bed->[vcf,idx]}
            .flatMap()
            .collect()
            .map{vcf_files->[[ id:(metadata.id?:"bam2vcf")],vcf_files.sort()]}
         )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

emit:
    versions
    vcf_chunks = concat
    vcf = BCFTOOLS_CONCAT.out.vcf
}

