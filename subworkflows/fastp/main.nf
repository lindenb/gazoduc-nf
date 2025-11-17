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
include { FASTQC as FASTQC_BEFORE             } from '../../modules/fastqc'
include { FASTQC as FASTQC_AFTER              } from '../../modules/fastqc'
include { FASTP  as APPLY_FASTP               } from '../../modules/fastp'
include { flatMapByIndex                      } from '../../modules/utils/functions'
workflow FASTP {
    take:
        metadata
        fastqs // [meta, [R1,R2]] or [meta,R0] or [meta,R1,R2] or [meta,[R0]]
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()
        fastqc = Channel.empty()

        


        fastqs = fastqs
                .map{
                    if(it.size()==3 && (it[1] instanceof Path) && (it[2]  instanceof Path)) {
                        return [it[0],[it[1],it[2]]];
                        }
                    if(it.size()==2 && (it[1] instanceof Path)) {
                        return [it[0],[it[1]]];
                        }
                    return it;
                    }
                .map {
                    if(it.size()!=2 || !(it[1] instanceof List)) {
                        throw new IllegalArgumentException("workflow:FASTP: Bad input ${it}.");
                        }
                    return it;
                    }

        fastqs_out = fastqs

         if(metadata.fastqc_before==null) {
                log.warn("FASTP: undefined metadata.fastqc_before")
                }
        if(metadata.fastqc_before==null || metadata.fastqc_before==true) {
                FASTQC_BEFORE( fastqs )
                versions = versions.mix(FASTQC_BEFORE.out.versions)
                fastqc = fastqc.mix(FASTQC_BEFORE.out.zip)
                }

        if(metadata.fastp_disabled==null || metadata.fastp_disabled==false ) {

           
            if(metadata.fastqc_after==null) {
                log.warn("FASTP: undefined metadata.fastqc_after")
                }
            
            
            APPLY_FASTP(  fastqs )
            versions = versions.mix(APPLY_FASTP.out.versions)
            multiqc = multiqc.mix(APPLY_FASTP.out.json)


            fastqs_out = APPLY_FASTP.out.fastqs

            if(metadata.fastqc_after==null || metadata.fastqc_after==true) {
                FASTQC_AFTER(fastqs_out )
                versions = versions.mix(FASTQC_AFTER.out.versions)
                fastqc = fastqc.mix(FASTQC_AFTER.out.zip)
                }
            }

            
        // split single and paired ends 
        trim_reads = fastqs_out.branch{v->
            paired_end : (v[1] instanceof List) && v[1].size()==2
            single_end:  (v[1] instanceof Path) || (v[1] instanceof List && v[1].size()==1)
            others: true
            }

        single_end = trim_reads.single_end
                .map{meta,fqs->[meta,(fqs instanceof List?fqs[0]:fqs)]}
                

        paired_end = trim_reads.paired_end
            .map{meta,fqs->[meta,fqs.sort()]}
            .map{meta,fqs->[meta,fqs[0],fqs[1]]}


        trim_reads.others.map{
                log.warn("FASTP:Illegal State $it .");
                throw new IllegalStateException("FASTP:Illegal State $it .");
                }
        

        multiqc = multiqc.mix(fastqc.flatMap{flatMapByIndex(it,1)})

    emit:
        versions
        multiqc
        single_end
        paired_end
        fastqs = fastqs_out //raw output of fastqs
}
