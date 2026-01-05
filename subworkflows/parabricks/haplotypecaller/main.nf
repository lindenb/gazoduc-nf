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
include {PB_HAPLOTYPECALLER as  PB_HAPCALLER   } from '../../../modules/parabricks/haplotypecaller'
include {GLNEXUS_GENOTYPE                      } from '../../../modules/glnexus/genotype'
include {BCFTOOLS_CONCAT                       } from '../../../modules/bcftools/concat3'

workflow PB_HAPLOTYPECALLER {
take:
    metadata
    fasta
    fai
    _dict
    cluster_beds // meta,bed
    glnexus_config
    bams // meta,bam,bai
main:
    versions = Channel.empty()
    multiqc = Channel.empty()

    PB_HAPCALLER(fasta,fai,bams)
    versions = versions.mix(PB_HAPCALLER.out.versions)

    ch1 = PB_HAPCALLER.out.gvcf
        .flatMap{_meta,vcf,tbi->[vcf,tbi]}//gvcf,tbi
        .collect()
        .map{files->[[id:(metadata.id?:"pb_hapcaller")],files.sort()]}
     
    GLNEXUS_GENOTYPE(
        cluster_beds,
        glnexus_config, //config
        ch1
        )
    versions = versions.mix(GLNEXUS_GENOTYPE.out.versions)

    BCFTOOLS_CONCAT(
        GLNEXUS_GENOTYPE.out.vcf
            .flatMap{_meta,vcf,tbi,_bed->[vcf,tbi]}//gvcf,tbi
            .flatMap()
            .collect()
            .map{files->[[id:(metadata.id?:"pb_hapcaller")],files.sort()]}
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

emit:
    versions
    multiqc
    vcf = BCFTOOLS_CONCAT.out.vcf
}
