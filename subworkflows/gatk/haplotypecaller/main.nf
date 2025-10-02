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

include { HAPLOTYPECALLER as HAPCALLER         }  from '../../../modules/gatk/hapcaller1'
include { BCFTOOLS_CONCAT                      }  from '../../../modules/bcftools/concat'
include { COMBINE_GENOTYPE_GVCFS               }  from '../combinegenotypegvcfs'
include { makeKey                              }  from '../../../modules/utils/functions.nf'
include {FIND_GVCF_BLOCKS                      }  from '../../..//modules/jvarkit/findgvcfblocks'


workflow HAPLOTYPECALLER {
take:
    meta
    fasta
    fai
    dict
    all_references// meta,fasta_files
    dbsnp // meta,vcf,tbi
    beds // meta,bed
    bams // meta,bam,bai
main:
    versions = Channel.empty()

    HAPCALLER(
        fasta,
        fai,
        dict,
        all_references,
        bams.combine(beds)
            .map{meta1,bam,bai,meta2,bed->[meta1,bam,bai,bed]}
        )
    versions  = versions.mix(HAPCALLER.out.versions)


    FIND_GVCF_BLOCKS(
        HAPCALLER.out.gvcf
            .map{meta,vcf,tbi,bed->[bed.toRealPath(),[vcf,tbi]]}
            .groupTuple()
            .map{bed,vcf_files->[[id:makeKey(bed)],vcf_files.flatten().sort(),bed]}
        )
    versions  = versions.mix(FIND_GVCF_BLOCKS.out.versions)

/*
    COMBINE_GENOTYPE_GVCFS(
            meta,
            fasta,
            fai,
            dict,
            HAPCALLER.out.gvcf
            )
    versions  = versions.mix(COMBINE_GENOTYPE_GVCFS.out.versions)


    BCFTOOLS_CONCAT(
        COMBINE_GENOTYPE_GVCFS.out.vcf
            .map{meta,vcf,tbi,bed->[vcf,tbi]}//gvcf,tbi
             .collect()
             .map{[[id:"gatkhapcaller"],it]},
        [[:],[]] //bed
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

emit:
    versions
    vcf = BCFTOOLS_CONCAT.out.vcf
*/
}
