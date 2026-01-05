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

include {BWA_MEM                         } from '../../../modules/bwa/mem'
include {SEQKIT_SPLIT                    } from '../../seqkit/split'
include {FASTP                           } from '../../../modules/fastp'
include {BQSR                            } from '../../gatk/bqsr'
include {MARK_DUPLICATES                 } from '../../markdup'
include {BAM_TO_CRAM                     } from '../../../modules/samtools/bam2cram'

workflow MAP_BWA {
	take:
		meta
		fasta
		fai
		dict
		BWADir
		bqsr_files
		bed
		fastqs
	main:
		versions = Channel.empty()
		out_bams = Channel.empty()
		multiqc  = Channel.empty()
		out_crams = Channel.empty()

		
		
		if(meta.with_fastp==null) {
			log.warn("undefined MAP_BWA:with_fastp")
			}

		
		fastqs = fastqs.map{
			if(it[1] instanceof List) {
				def L=it[1].sort();
				if(L.size()==2) return [it[0], L[0], L[1]];
				if(L.size()==1) return [it[0], L[0], null];
				throw new IllegalArgumentException("MAP_BWA : bad fastq input : ${it}");
				}
			return it;
			}

		if(meta.with_fastp==null || meta.with_fastp==true) {
			FASTP(
				fastqs.map{met,R1,R2->{
					if(R2) return [met,[R1,R2]];
					return [met,[R1]]
					}}
				)
			versions = versions.mix(FASTP.out.versions)
			

			fastqs = FASTP.out.fastqs.map{meta2,fqs->{
					if(fqs.size()==1) return [meta2,fqs[0],[]];
					if(fqs.size()!=2) throw new IllegalArgumentException("Boum after FASTP ${meta2} ${fqs}"); 
					def L1 = fqs.sort();
					return [meta2,L1[0],L1[1] ];
					}
				}
			}
		
		
		if(meta.with_seqkit_split==null || meta.with_seqkit_split==true) {
			SEQKIT_SPLIT(
				meta,
				fastqs
				)
			versions = versions.mix(SEQKIT_SPLIT.out.versions)
			fastqs = SEQKIT_SPLIT.out.fastqs
			}
		

		
		BWA_MEM(fasta,fai,BWADir,bed,fastqs)
		versions = versions.mix(BWA_MEM.out.versions)
		
		
		MARK_DUPLICATES(
			meta.plus(markdup_enabled : meta.with_markdup),
			fasta,
			fai,
			dict,
			BWA_MEM.out.bam
			)
		versions = versions.mix(MARK_DUPLICATES.out.versions)
		multiqc = multiqc.mix(MARK_DUPLICATES.out.multiqc)
		out_bams = MARK_DUPLICATES.out.bams


		if(meta.with_bqsr==null) {
			log.warn("undefined MAP_BWA:with_bqsr");
			}
		if(meta.with_bqsr==null || meta.with_bqsr==true) {
			BQSR(
				meta,
				fasta,
				fai,
				dict,
				bqsr_files,
				out_bams
				)
			versions = versions.mix(BQSR.out.versions)
			out_bams = BQSR.out.bams
			}
			

		if(meta.with_cram==null) {
			log.warn("MAP_BWA: undefined meta.with_cram")
			}


		if(meta.with_cram==null || meta.with_cram==true) {
			BAM_TO_CRAM(fasta,fai,out_bams.map{meta,bam,bai->[meta,bam]});
			versions = versions.mix(BAM_TO_CRAM.out.versions)
			out_crams = out_crams.mix(BAM_TO_CRAM.out.cram)
			}
		
	emit:
		versions
		bams = out_bams
		crams = out_crams
		multiqc
	}
