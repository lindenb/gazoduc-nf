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
include { GATK_BAM2VCF   } from '../../../modules/gatk/bam2vcf'
include { makeKey        } from '../../../modules/utils/functions'
include { BED_SPLITX     } from '../../../modules/bed/splitx'



workflow BAM2VCF_DIVIDE_AND_CONQUER {
take:
    workflow_meta
    level
    fasta
    fai
    dict
    dbsnp
    pedigree
    references
    bed // [meta,bed]
    bams // [meta, [bams and bai] ]
main:
    versions = Channel.empty()
    bed_todo = bed.map{meta,bed->[meta.plus("level":level),bed]}

    GATK_BAM2VCF(
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        bams.combine(bed_todo)
            .map{_meta1,bam_files,meta2,bed->[
                meta2/* bed.meta */,
                bam_files/*bam and bai */,
                bed/*bed */
                ]}
        )
    versions = versions.mix(GATK_BAM2VCF.out.versions)

    bed_todo 
        .mix(GATK_BAM2VCF.out.vcf.map{meta,_vcf,_tbi,bed->[meta,bed]})
        .groupTuple()
        .branch{_meta,beds->
            success: beds.size()==2
            failure: beds.size()==1
            other: true
            }.set{branch2}
            
    branch2.other.map{throw new IllegalArgumentException(" BAM2VCF_DIVIDE_AND_CONQUER ${it}");}

    BED_SPLITX(
        branch2.failure
        )
    versions = versions.mix(BED_SPLITX.out.versions)

    bed_todo = BED_SPLITX.out.beds
        .map{_meta,beds->beds}
        .map{beds->(beds instanceof List?beds:[beds])}
        .flatMap()
        .map{bed->[[id:makeKey(bed)],bed]}
       
    vcf_out = GATK_BAM2VCF.out.vcf

emit:
    versions
    bed = bed_todo
    vcf = vcf_out
}


