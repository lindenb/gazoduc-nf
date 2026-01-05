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
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER1 } from './sub.nf'
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER2 } from './sub.nf'
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER3 } from './sub.nf'
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER4 } from './sub.nf'
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER5 } from './sub.nf'
include {BAM2VCF_DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER6 } from './sub.nf'
include {BCFTOOLS_CONCAT                                    } from '../../../modules/bcftools/concat3'



/** Divide an conquer BAM2vcf */
workflow GATK_BAM2VCF_DNC {
take:
    workflow_meta
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
    concat = Channel.empty()

    bams_ch = bams
        .map{meta,bam,bai->[bam,bai]}
        .flatMap()
        .collect()
        .map{[[id: (workflow_meta.id?:"bam2vcf.dnc")],it.sort()]}
    

    DIVIDE_AND_CONQUER1(
        workflow_meta,
        1,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER1.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER1.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER1.out.versions)


    DIVIDE_AND_CONQUER2(
        workflow_meta,
        2,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER2.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER2.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER2.out.versions)

    DIVIDE_AND_CONQUER3(
        workflow_meta,
        3,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER3.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER3.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER3.out.versions)

    DIVIDE_AND_CONQUER4(
        workflow_meta,
        4,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER4.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER4.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER4.out.versions)

    DIVIDE_AND_CONQUER5(
        workflow_meta,
        5,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER5.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER5.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER5.out.versions)

    DIVIDE_AND_CONQUER6(
        workflow_meta,
        6,
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER6.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER6.out.vcf)
    versions = versions.mix(DIVIDE_AND_CONQUER6.out.versions)

 

      
    BCFTOOLS_CONCAT(
        concat.map{_meta,vcf,tbi,_bed->[vcf,tbi]}
            .flatMap()
            .collect()
            .map{[[id:(workflow_meta.id?:"gatk.dnc")],it.sort()]}
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)
    out_vcf = BCFTOOLS_CONCAT.out.vcf
emit:
    versions
    vcf_chunks = concat /* meta,vcf,tbi,bed */
    vcf = out_vcf
}
