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

include {BWA_MEM                         } from '../../../modules/bwa/mem'
include {SEQKIT_SPLIT                    } from '../../seqkit/split'
include {FASTP                           } from '../../../modules/fastp'
include {MARK_DUPLICATES                 } from '../../../modules/gatk/markduplicates'
include {BQSR                            } from '../../gatk/bqsr'

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
		
		fastqs.view()
		
		
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
		
		fastqs.view{"after fastp ${it}"}
	
		if(meta.with_seqkit_split==null || meta.with_seqkit_split==true) {
			SEQKIT_SPLIT(
				meta,
				fastqs
				)
			versions = versions.mix(SEQKIT_SPLIT.out.versions)
			fastqs = SEQKIT_SPLIT.out.fastqs
			}
		
		fastqs.view{"after seqkit ${it}"}

		BWA_MEM(fasta,fai,BWADir,bed,fastqs)
		versions = versions.mix(BWA_MEM.out.versions)
	
		if(meta.with_markdup==null || meta.with_markdup==true) {
			if(meta.markdup==null || meta.markdup.equals("markduplicates")) {
				MARK_DUPLICATES(BWA_MEM.out.bam
					.map{meta,bam,bai->[meta.id,meta,[bam,bai]]}
					.groupTuple().view()
					.map{id,metas,bam_files->[metas[0],bam_files.flatten().sort()]}
					)
				versions = versions.mix(MARK_DUPLICATES.out.versions)
				out_bams = MARK_DUPLICATES.out.bam
				}
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
		
	emit:
		versions
		bam = out_bams
	}
