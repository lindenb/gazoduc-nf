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

include {MANTA_SINGLE as MANTA     } from '../../../modules/manta/single'
include {MANTA_CANDIDATESV         } from '../../../modules/manta/candidateSV'
include {MANTA_CONVERT_INVERSION   } from '../../../modules//manta/convert.inversion'
include {TRUVARI_COLLAPSE          } from '../../../modules/truvari/collapse'

workflow MANTA_SINGLE {
take:
	meta
        fasta
        fai
        dict
        bams
	optional_bed
main:
        versions = Channel.empty()

        
	
        MANTA(
                fasta,
                fai,
                dict,
                optional_bed,
                bams
                )

        versions = versions.mix(MANTA.out.versions)
        
        MANTA_CONVERT_INVERSION(
				fasta,
				fai,
				dict,
                MANTA.out.diploidSV,
                )
        versions = versions.mix(MANTA_CONVERT_INVERSION.out.versions)




        ch1 = MANTA_CONVERT_INVERSION.out.vcf
                .map{[it[1],it[2]]}
                .flatten()
                .collect()
                .filter{it.size()>2} // NO need to run truvari if there is only one vcf and one tbi
                .map{[[id:"manta"],it]}
                


        MANTA_CANDIDATESV( MANTA.out.candidateSV)
        versions = versions.mix(MANTA_CANDIDATESV.out.versions)
        

	if(meta.with_truvari==null) {
		log.warn("MANTA: meta.with_truvari undefined")
		}        

	if(meta.with_truvari==null || meta.with_truvari==true) {
                fasta2  = bams
                        .count()
                        .filter{it>1}
                        .combine(fasta)
                        .map{_count,meta,fasta->[meta,fasta]}
                        .first()
                        
	        TRUVARI_COLLAPSE(
			fasta2,
			fai,
			dict,
			ch1
			)
        	versions = versions.mix(TRUVARI_COLLAPSE.out.versions)
		}
emit:
        versions
}
