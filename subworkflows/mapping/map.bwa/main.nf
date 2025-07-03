/*

Copyright (c) 2024 Pierre Lindenbaum

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

include {BWA_MEM_01} from '../../modules/bwa/bwa.mem.01.nf'
include {SAMBAMBA_MARKDUP_01} from '../../modules/sambamba/sambamba.markdup.01.nf'
include {SAMTOOLS_MERGE_01} from '../../modules/samtools/samtools.merge.01.nf'
include {SAMTOOLS_BAM_TO_CRAM_01} from '../../modules/samtools/samtools.bam2cram.01.nf'
include {MERGE_VERSION as MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {moduleLoad;isBlank;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {SEQTK_SPLITFASTQ} from '../fastq/split.fastq.nf'
include {ORA_TO_FASTQ}  from '../ora/ora2fastq.nf'

workflow MAP_BWA_01 {
	take:
		meta
		fasta
		fai
		dict
		fastqs
	main:
		versions = Channel.empty()
		def max_size = meta.max_size?:1_000_000_000L;
	
		choice = fastqs.branch {
			big: it[1].size() >= max_size || (itsize()>2 && it[2].size() >= max_size)
			small: true
		}

		fastq_ch2 =  SEQTK_SPLITFASTQ(choice.big)
		versions = version_ch.mix(fastq_ch2.versions)


		R1 = SEQTK_SPLITFASTQ.out.R1
			.flatMap{T->T[1].map{[T[0],it]}}
		R2 = SEQTK_SPLITFASTQ.out.R2
			flatMap{T->T[1].map{[T[0],it]}}

		r1r2_ch = R1.join(R2)

	
		bam1_ch = BWA_MEM_01(fasta,fai,dict, r1r2_ch.mix(choice.small))
		versions = version_ch.mix(bam1_ch.versions)
	

		merge_ch = SAMTOOLS_MERGE_01([:], genomeId, bam1_ch.output.map{T->[T[0].sample,T[1]] }.groupTuple())
		version_ch = version_ch.mix(merge_ch.version)

		markdup_ch = SAMBAMBA_MARKDUP_01([:],merge_ch.bam)
		version_ch = version_ch.mix(markdup_ch.version)

		if( parseBoolean(params.bwa.with_bqsr) ) {
			
		
			cram_ch = SAMTOOLS_BAM_TO_CRAM_01([:],applybqsr_ch.bam.map{T->["sample":T[0],"bam":T[1],"genomeId":genomeId]})
			version_ch = version_ch.mix(cram_ch.version)
			}
		else
			{
			cram_ch = SAMTOOLS_BAM_TO_CRAM_01([:], markdup_ch.bam.map{T->["sample":T[0],"bam":T[1],"genomeId":genomeId]})
			version_ch = version_ch.mix(cram_ch.version)		
			}

		version_ch = MERGE_VERSION("bwa", version_ch.collect())
	emit:
		version= version_ch.version
		bams = cram_ch.output
	}
