/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {SCATTER_TO_BED} from '../../subworkflows/picard/picard.scatter2bed.nf'
include {SAMTOOLS_FASTQ_01} from '../../modules/samtools/samtools.collate.fastq.01.nf'
include {BWA_MEM_01} from '../../modules/bwa/bwa.mem.01.nf'


workflow REMAP_BWA_01 {
	take:
		meta
		reference_in
		reference_out
		bams
	main:
		version_ch = Channel.empty()

		acgt_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"100000"],reference_in)
		version_ch = version_ch.mix(acgt_ch.version)


		all_samples_ch = SAMTOOLS_SAMPLES01([:],reference_in,bams)
		version_ch = version_ch.mix(all_samples_ch.version)

		each_sn_bam =  all_samples_ch.output.splitCsv(header:false,sep:'\t')

		intervals_ch = acgt_ch.bed.splitCsv(header:false,sep:'\t').
				map{T->T[0]+":"+((T[1] as int)+1)+"-"+T[2]}.
				mix(Channel.of("*"))


		remap_ch = REMAP_ONE(meta, reference_in, reference_out, each_sn_bam.combine(intervals_ch))
		version_ch = version_ch.mix(remap_ch.version)

	emit:
		version = version_ch
	}

workflow REMAP_ONE {
take:
	meta
	reference_in
	reference_out
	sample_bam_interval
main:
	rows_ch = sample_bam_interval.map{T->[
		"sample": T[0],
		"bam": T[1],
		"reference": reference_in,
		"interval": T[2]
		]}

	fastq_ch = SAMTOOLS_FASTQ_01(meta,rows_ch)

	map1_ch = fastq_ch.output.splitCsv(header:true,sep:'\t')

	umap1_ch = map1_ch.flatMap{T->[
		[[T.sample,"R1"],T.unpairedR1],
		[[T.sample,"R2"],T.unpairedR2]			
		]}.groupTuple()
	unmap2_ch = SORT_UNPAIRED_FASTQ(meta, umap1_ch)
	
	unmap3_R1_ch = unmap2_ch.filter{it[1].equals("R1")}.map{T->[T[0],T[2]]}
	unmap3_R2_ch = unmap2_ch.filter{it[1].equals("R2")}.map{T->[T[0],T[2]]}
	unmap4_ch = JOIN_UNPAIRED(meta,unmap3_R1_ch.join(unmap3_R2_ch))


	r1r2_ch= map1_ch.map{T->[
		"sample":T.sample,
		"R1":T.R1,
		"R2":T.R2,
		"reference":reference_out
		]}
	
	bam2_ch = unmap4_ch.output.map{T->[
                "sample":T[0],  
                "R1":T[1],
		"interleaved":true,
                "reference":reference_out
                ]}
	

	
	bam1_ch = BWA_MEM_01(meta,r1r2_ch.mix(bam2_ch))

	merge_ch = MERGE_BAMS(meta, bam1_ch.bam.groupTuple())

emit:
	version= fastq_ch.version
}


process SORT_UNPAIRED_FASTQ {
tag "${key[0]} ${key[1]} N=${L.size()}"
afterScript "rm -f jeter.tsv.gz"
memory "3g"
input:
	val(meta)
	tuple val(key),val(L)
output:
	tuple val("${key[0]}"),val("${key[1]}"),path("${key[0]}.${key[1]}.tsv.gz"),emit:output
script:
"""
gunzip -c ${L.join(" ")} |\
	LC_ALL=C sort -T . -S ${task.memory.kilo} -t '\t' -k1,1 |\
	gzip --best > jeter.tsv.gz

mv jeter.tsv.gz "${key[0]}.${key[1]}.tsv.gz"
"""
}

process JOIN_UNPAIRED {
tag "${sample} ${R1} ${R2}"
afterScript "rm -f jeter.fq.gz jeter2.fq.gz"
memory "3g"
input:
	val(meta)
	tuple val(sample),val(R1),val(R2)
output:
	tuple val(sample),path("${sample}.R1R2.fq.gz"),emit:output
	tuple val(sample),path("${sample}.R0.fq.gz"),emit:single
script:
"""


LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,1.4,2.1,2.2,2.3,2.4' \
	<(gunzip -c "${R1}")  \
	<(gunzip -c "${R2}") |\
	tr "\t" "\\n" |\
	gzip --best > jeter.fq.gz

LC_ALL=C join -t '\t' -1 1 -2 1 -v 1 -v2 \
	<(gunzip -c "${R1}")  \
	<(gunzip -c "${R2}") |\
	tr "\t" "\\n" |\
	gzip --best > jeter2.fq.gz


mv jeter.fq.gz  "${sample}.R1R2.fq.gz"
mv jeter2.fq.gz "${sample}.R0.fq.gz"
"""
}

process MERGE_BAMS {
tag "${sample}"
tag "N=${L.size()}"
input:
	val(meta)
	tuple val(sample),val(L)
output:
	tuple val(sample),path("${sample}.merged.bam"),emit:bam
script:
"""
module load samtools

samtools merge -o "${sample}.merged.bam" ${L.join(" ")}
"""
}
