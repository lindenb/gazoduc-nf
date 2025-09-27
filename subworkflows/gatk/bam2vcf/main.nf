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
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER1 } from './sub.nf'
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER2 } from './sub.nf'
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER3 } from './sub.nf'
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER4 } from './sub.nf'
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER5 } from './sub.nf'
include {DIVIDE_AND_CONQUER as  DIVIDE_AND_CONQUER6 } from './sub.nf'


workflow GATK_BAM2VCF {
take:
    meta
    fasta
    fai
    dict
    dbsnp
    references //[meta, [ref files fa fai dict...]] all known reference
    beds // [meta,bed]
    bams // [meta,bam,bai]
main:
    versions = Channel.empty()
    concat = Channel.empty()

    bams_ch = bams
        .map{[it[1],it[2]]}
        .flatMap()
        .collect()
        .map{[[id:"hapcaller"],it.sort()]}
    

    DIVIDE_AND_CONQUER1(
        meta,
        1,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER1.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER1.out.vcf)


    DIVIDE_AND_CONQUER2(
        meta,
        2,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER2.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER2.out.vcf)

    DIVIDE_AND_CONQUER3(
        meta,
        3,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER3.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER3.out.vcf)


    DIVIDE_AND_CONQUER4(
        meta,
        4,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER4.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER4.out.vcf)


    DIVIDE_AND_CONQUER5(
        meta,
        5,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER5.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER5.out.vcf)

    DIVIDE_AND_CONQUER6(
        meta,
        6,
        fasta,
        fai,
        dict,
        dbsnp,
        references,
        beds,
        bams_ch
        )
    beds = DIVIDE_AND_CONQUER6.out.bed
    concat = concat.mix(DIVIDE_AND_CONQUER6.out.vcf)



emit:
    versions
}