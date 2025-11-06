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
include { BAM_TO_FASTQ as BAM2FASTQ } from '../../modules/samtools/bam2fastq'
include { SAM_TO_FASTQ              } from '../../modules/gatk/sam2fastq'
include { PB_BAM2FQ                 } from '../../modules/parabricks/bam2fq'
include { isEmptyGz                 } from '../../modules/utils/functions.nf'



workflow BAM_TO_FASTQ {
take:
    workflow_metadata // bam2fastq_method=(samtools|gatk|parabricks)
    bams // meta,bam,bai,fasta,fai,dict,opt_bed
main:
    versions = Channel.empty()
    paired = Channel.empty()
    single =  Channel.empty()
    
    /* convert to meta, bam, bai,fasta,fai,dict,opt_bed if needed */
    ch1= bams.map{
        if(it.size()==7) {
            return it;
            }
        def L=[];
        L.addAll(it);
        while(L.size()<7) L.add([]);
        return L;
        }
    if(workflow_metadata.bam2fastq_method==null) {
        log.warn("BAM_TO_FASTQ: undefined bam2fastq_method. using samtools.");
        }
    /* ============================================================================================
     *
     * SAMTOOLS
     *
     */
    if(workflow_metadata.bam2fastq_method==null || workflow_metadata.bam2fastq_method.equalsIgnoreCase("samtools")) {
        BAM2FASTQ(ch1.map{meta,bam,bai,fasta,fai,_dict,bed->[meta,bam,bai,fasta,fai,bed]})
        versions = versions.mix(BAM2FASTQ.out.versions)
        paired = paired.mix(BAM2FASTQ.out.fastq.map{meta,R1,R2,R0,rS->[meta,R1,R2]})
        single = single
            .mix(BAM2FASTQ.out.fastq.map{meta,R1,R2,R0,RS->[meta,R0]})
            .mix(BAM2FASTQ.out.fastq.map{meta,R1,R2,R0,RS->[meta,RS]})
        }
    /* ============================================================================================
     *
     * GATK/PICARD
     *
     */
    else if(workflow_metadata.bam2fastq_method.equalsIgnoreCase("gatk") || workflow_metadata.bam2fastq_method.equalsIgnoreCase("picard")) {
        // BED is not supported by picard or gatk
        ch1.map{meta,bam,bai,fasta,fai,dict,bed->[meta,bed]}
            .filter{meta,bed->!(bed instanceof List && bed.isEmpty())}
            .map{meta,bed->throw new IllegalArgumentException("BAM_TO_FASTQ:BED not supported gatk and ${meta}")}

        ch2 = ch1.multiMap{meta,bam,bai,fasta,fai,dict,_bed->
            fasta: [meta,fasta]
            fai: [meta,fai]
            dict: [meta,dict]
            bam: [meta,bam,bai]
            }
        SAM_TO_FASTQ(
            ch2.fasta,
            ch2.fai,
            ch2.dict,
            ch2.bam
            )
        versions = versions.mix(SAM_TO_FASTQ.out.versions)
        paired = paired.mix(SAM_TO_FASTQ.out.paired_end.map{meta,R1,R2,R0->[meta,R1,R2]})
        single = single
            .mix(SAM_TO_FASTQ.out.single_end)
            .mix(SAM_TO_FASTQ.out.paired_end.map{meta,R1,R2,R0->[meta,R0]})
        }
    /* ============================================================================================
     *
     * PARABRICKS
     *
     */
    else if(workflow_metadata.bam2fastq_method.equalsIgnoreCase("parabricks")) {
        // BED is not supported by picard or gatk
        ch1.map{meta,bam,bai,fasta,fai,dict,bed->[meta,bed]}
            .filter{meta,bed->!(bed instanceof List && bed.isEmpty())}
            .map{meta,bed->throw new IllegalArgumentException("BAM_TO_FASTQ:BED not supported parabrick and ${meta}")}

        ch2 = ch1.multiMap{meta,bam,bai,fasta,fai,dict,_bed->
            fasta: [meta,fasta]
            fai: [meta,fai]
            dict: [meta,dict]
            bam: [meta,bam,bai]
            }
        
        PB_BAM2FQ( 
            ch1.map{meta,bam,_bai,fasta,fai,_dict,_bed->[meta,bam,fasta,fai]}
            )
        versions = versions.mix(PB_BAM2FQ.out.versions)
        paired = paired.mix(
            PB_BAM2FQ.out.fastqs
                .filter{meta,fastqs->[meta,(fastqs instanceof List?fastqs:[fastqs])]} //may be single end
                .filter{meta,fastqs->fastqs.size()==2} //check paired end
                .map{meta,fastqs->[meta,fastqs.sort()]} // check R1 before R2
                .map{meta,fastqs->[meta,fastqs[0],fastqs[1]]}
            )
        single = single.mix(
            PB_BAM2FQ.out.fastqs
                .filter{meta,fastqs->[meta,(fastqs instanceof List?fastqs:[fastqs])]} //may be single end
                .filter{meta,fastqs->fastqs.size()==1} //check single end
                .map{meta,fastqs->[meta,fastqs[0]]}
            )
        }
    else
        {
        throw new IllegalArgumentException("undefined BAM_TO_FASTQ.bam2fastq_method= ${workflow_metadata.bam2fastq_method}")
        }

    paired = paired.filter{meta,R1,R2->!(isEmptyGz(R1) && isEmptyGz(R2))}
    single = single.filter{meta,R1->!(isEmptyGz(R1))}
emit:
    versions
    paired_end = paired
    single_end = single
}
