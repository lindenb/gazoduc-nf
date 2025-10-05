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
			fastqs = FASTP.out.fastqs
			}
		
		
		
		if(meta.with_seqkit_split==null || meta.with_seqkit_split==true) {
			SEQKIT_SPLIT(meta,fastqs.map{meta,fqs->{
				if(fqs.size()==1) return [meta,fqs[0],[]];
				if(fqs.size()!=2) throw new IllegalArgumentException("Boum after FASTP"); 
				def L1 = fqs.sort();
				return [meta,L1[0],L1[1]];
				}})
			versions = versions.mix(SEQKIT_SPLIT.out.versions)
			fastqs = SEQKIT_SPLIT.out.fastqs
			}
		
		
		BWA_MEM(fasta,fai,BWADir,bed, 
			fastqs.map{meta,fqs->{
				if(fqs.size()==1) return [meta,fqs[0],[]];
				if(fqs.size()!=2) throw new IllegalArgumentException("Boum after FASTP"); 
				def L1 = fqs.sort();
				return [meta,fqs[0],fqs[1]];
				}})
		versions = versions.mix(BWA_MEM.out.versions)
	
		
	
		if(meta.markdup==null || meta.markdup.equals("markduplicates")) {
				MARK_DUPLICATES(BWA_MEM.out.bam
					.map{meta,bam,bai->[meta.id,meta,[bam,bai]]}
					.groupTuple().view()
					.map{id,metas,bam_files->[metas[0],bam_files.flatten.sort()]}
				)
			versions = versions.mix(MARK_DUPLICATES.out.versions)
			out_bams = MARK_DUPLICATES.out.bam
			}

		
	emit:
		versions
		bam = out_bams
	}
