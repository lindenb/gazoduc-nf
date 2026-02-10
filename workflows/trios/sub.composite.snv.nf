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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



*/
include { BCFTOOLS_CONCAT                          } from '../../modules/bcftools/concat3/main.nf'
include { HET_COMPOSITE                            } from '../../modules/jvarkit/hetcomposite/main.nf'

workflow WORKFLOW_COMPOSITE_SNV {
take:
    metadata
    fasta
    fai
    dict
    gff3
    gtf
    triosbams_ch
    pedigree
    vcf
main:
    versions = Channel.empty()

     BCFTOOLS_CONCAT(
        vcf.map{_meta,vcf,idx->[vcf,idx]}
            .flatMap()
            .collect()
            .map{files->[[id:"hetcomposite"],files.sort()]}
        )
    versions =  versions.mix(BCFTOOLS_CONCAT.out.versions)
    vcf = BCFTOOLS_CONCAT.out.vcf

    HET_COMPOSITE(
        fasta,
        fai,
        dict,
        pedigree,
        vcf
        )
    versions =  versions.mix(HET_COMPOSITE.out.versions)
emit:
    versions
}