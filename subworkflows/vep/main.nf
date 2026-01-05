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
include {VEP as VEP_GRCH38 } from './grch38.nf'
include {isBlank           } from '../../modules/utils/functions.nf'

workflow VEP {
    take:
        metadata
        fasta
        fai
        dict
        vcf
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()
        out_vcf = Channel.empty()
        ch1 = fasta.branch{meta,_fasta->
            fasta_hg38: !isBlank(meta.ucsc_name) && meta.ucsc_name=="hg38"
            other: true
            }


        VEP_GRCH38(
            metadata,
            fasta_hg38,
            fai,
            dict,
            vcf
            )
        out_vcf = out_vcf.mix(VEP_GRCH38.out.vcf)
        // out_vcf = out_vcf.mix(vcf)



    emit:
        versions
        multiqc
        vcf
}