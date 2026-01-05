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
include { XHUNTER_DENOVO_PROFILE       } from '../../../modules/expansion.hunter/denovo.profile'
include { XHUNTER_DENOVO_MERGE         } from '../../../modules/expansion.hunter/denovo.merge'
include { XHUNTER_DENOVO_CASE_CONTROL  } from '../../../modules/expansion.hunter/denovo.case.control'
include { isBlank                      } from '../../../modules/utils/functions.nf'

workflow EXPANSION_HUNTER_DE_NOVO {
	take:
		workflow_metadata
		fasta
        fai
        dict
        bams
	main:
        versions = Channel.empty()

        ch1 = bams.branch{
            bam3: it.size()==3
            bam5: it.size()==5
            other: true
        }

        ch1.other.map{throw new IllegalArgumentException("EXPANSION_HUNTER_DE_NOVO bad input");}

        ch2 = ch1.bam5.mix(
            ch1.bam3.combine(fasta)
                .combine(fai)
                .map{meta1,bam,bai,_meta2,fasta,_meta3,fai->[meta1,bam,bai,fasta,fai]}
            ).multiMap{meta,bam,bai,fasta,fai->
                fasta: [meta,fasta]
                fai: [meta,fai]
                bam: [meta,bam,bai]
            }


        XHUNTER_DENOVO_PROFILE(
            ch2.fasta,
            ch2.fai,
            ch2.bam
            )
		versions = versions.mix(XHUNTER_DENOVO_PROFILE.out.versions)

        json = XHUNTER_DENOVO_PROFILE.out.json
            .filter{meta,_json->!isBlank(meta.status)}
        

        manifest_ch = json
                .map{meta,json->"${meta.id}\t${meta.status}\tJSON/${json.name}"}
                .collectFile(name:"manifest.tsv",newLine:true)
                .map{[[id:workflow_metadata.id],it]}


		XHUNTER_DENOVO_MERGE(
            fasta,
            fai,
            json
                .map{_meta,json->json}
                .collect()
                .map{[[id:workflow_metadata.id],it]},
            manifest_ch
            )

		versions = versions.mix(XHUNTER_DENOVO_MERGE.out.versions)

		methods = Channel.of("motif","locus")

		XHUNTER_DENOVO_CASE_CONTROL(
            fasta,
            fai,
            XHUNTER_DENOVO_MERGE.out.json,
            manifest_ch
                .combine(methods)
                .map{meta,manifest,meth->[meta.plus(method:meth),manifest]}
            )
        versions = versions.mix(XHUNTER_DENOVO_CASE_CONTROL.out.versions)
	emit:
		versions
        tsv = XHUNTER_DENOVO_CASE_CONTROL.out.tsv
	}






