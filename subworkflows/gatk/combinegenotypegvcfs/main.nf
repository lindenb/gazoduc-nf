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
include {COMBINE_GVCFS     } from '../combinegvcfs'
include {GENOTYPEGVCFS     } from '../../../modules/gatk/genotypegvcfs'

workflow COMBINE_GENOTYPE_GVCFS {
    take:
        meta
        fasta
        fai
        dict
        dbsnp // meta, vcf,tbi
        gvcfs_bed // [meta, gvcfgz, gvcfgz_tbi, bed]
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()

        COMBINE_GVCFS(meta,fasta,fai,dict,gvcfs_bed)
        versions = versions.mix(COMBINE_GVCFS.out.versions)
        multiqc = multiqc.mix(COMBINE_GVCFS.out.multiqc)

        GENOTYPEGVCFS(fasta,fai,dict,dbsnp, COMBINE_GVCFS.out.gvcf)
        versions = versions.mix(GENOTYPEGVCFS.out.versions)

        vcf = GENOTYPEGVCFS.out.vcf

    emit:
        versions
        multiqc
        vcf

}
