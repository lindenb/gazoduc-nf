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
include {GATK_BAM2VCF   } from '../../../modules/gatk/bam2vcf'
include {makeKey        } from '../../../modules/utils/functions'



workflow BAM2VCF_DIVIDE_AND_CONQUER {
take:
    meta
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
    ok_ch = Channel.empty()
    bed_todo = bed.map{meta,bed->[meta.plus("level":level),bed]}

   

    GATK_BAM2VCF(
        fasta,
        fai,
        dict,
        dbsnp,
        pedigree,
        references,
        bams.combine(bed_todo)
            .map{meta1,bam_files,meta2,bed->[
                meta2/* bed.meta */,
                bam_files/*bam and bai */,
                bed/*bed */
                ]}
        )
    versions = versions.mix(GATK_BAM2VCF.out.versions)
    

    bed_todo 
        .mix(GATK_BAM2VCF.out.vcf.map{meta,vcf,tbi,bed->[meta,bed]})
        .groupTuple()
        .branch{meta,beds->
            success: beds.size()==2
            failure: beds.size()==1
            other: true
            }.set{branch2}
    branch2.other.map{throw new IllegalArgumentException("${it}");}


    SPLITBED(
        fasta,
        fai,
        dict,
        level,
        branch2.failure
        )
    
    bed_todo = SPLITBED.out.bed
        .map{meta,beds->beds}
        .map{beds instanceof List?beds:[beds]}
        .flatMap()
        .map{bed->[[id:makeKey(bed)],bed]}
       
    vcf_out =  GATK_BAM2VCF.out.vcf

   vcf_out.view()
emit:
    versions
    bed = bed_todo
    vcf = vcf_out
}


process SPLITBED {
tag "${bed.name} Level ${level}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    val(level)
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("BEDS/*"),optional:true,emit:bed
    path("versions.yml"),emit:versions
script:
    def njobs  = task.ext.njobs?:"10"
"""
python3 ${moduleDir}/../../../src/python/split_bed.py" ${njobs} '${bed}'

touch versions.yml
"""
}
