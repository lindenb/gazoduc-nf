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
include {SAMTOOLS_MERGE                  } from '../../../modules/samtools/merge'
include {BQSR                            } from '../../gatk/bqsr'
include {SAMBAMBA_MARKDUP                } from '../../../modules/sambamba/markdup'
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
		multiqc = Channel.empty()

		
		
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
		
		
	
		if(meta.with_markdup==null || meta.with_markdup==true) {
			if(meta.markdup_method==null || meta.markdup_method.equals("markduplicates")) {
				MARK_DUPLICATES(BWA_MEM.out.bam
					.map{meta,bam,bai->[meta.id,meta,[bam,bai]]}
					.groupTuple()
					.map{id,metas,bam_files->[metas[0],bam_files.flatten().sort()]}
					)
				versions = versions.mix(MARK_DUPLICATES.out.versions)
				out_bams = MARK_DUPLICATES.out.bam
				multiqc = multiqc.mix(MARK_DUPLICATES.out.metrics)
				}
			else if(meta.markdup_method.equals("sambamba")) {
			
			
				ch1 = BWA_MEM.out.bam
					.map{meta,bam,bai->[meta.id,meta,[bam,bai]]}
					.groupTuple()
					.map{id,metas,bam_files->[metas[0],bam_files]}
					.branch{
						needmerge: it[1].size() > 1
						one: true
						}
				
				
				SAMTOOLS_MERGE(
					fasta,
					fai,
					ch1.needmerge.map{meta,bam_files->[meta,bam_files.flatten().sort()]}
					)
				versions = versions.mix(SAMTOOLS_MERGE.out.versions)
				
				
				SAMBAMBA_MARKDUP(					
					SAMTOOLS_MERGE.out.bam.mix(ch1.one.map{meta,bam_files->[meta,bam_files[0]]})
					)
				versions = versions.mix(SAMBAMBA_MARKDUP.out.versions)
				out_bams = SAMBAMBA_MARKDUP.out.bams
				}
			else
				{
				throw new IllegalArgumentException("Boum MAP_BWA undefined meta.markdup_method: ${meta.markdup_method}");
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
			
		out_crams  =Channel.empty()
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
