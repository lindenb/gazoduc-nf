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
include { XHUNTER_APPLY                } from '../../modules/expansion.hunter/apply'
include { BCFTOOLS_MERGE               } from '../../modules/bcftools/merge'
include { isBlank                      } from '../../modules/utils/functions.nf'

workflow EXPANSION_HUNTER {
	take:
		workflow_metadata
		fasta
        fai
        dict
        catalog
        bams
	main:
        versions = Channel.empty()



        XHUNTER_APPLY(
            fasta,
            fai,
            catalog,
            bams
            )
		versions = versions.mix(XHUNTER_APPLY.out.versions)

        BCFTOOLS_MERGE(
            XHUNTER_APPLY.out.vcf
                .map{meta,vcf,tbi->[vcf,tbi]}
                .flatMap()
                .collect()
                .map{[
                    [id:workflow_metadata.id],
                    it.sort(),
                    []
                    ]}
            )
        versions = versions.mix(BCFTOOLS_MERGE.out.versions)
        
	emit:
		versions
        vcf = BCFTOOLS_MERGE.out.vcf
	}






