/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {MERGE_VERSION as MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GATK4_BASE_RECALIBRATOR_01} from '../../modules/gatk/gatk4.base.recalibration.01.nf'
include {GATK4_GATHER_BQSR_01} from '../../modules/gatk/gatk4.gather.bqsr.01.nf'
include {GATK4_APPLY_BQSR_01} from '../../modules/gatk/gatk4.apply.bqsr.01.nf'
include {moduleLoad;isBlank;parseBoolean} from '../../modules/utils/functions.nf'

workflow MAP_BWA_01 {
	take:
		meta
		reference
		fastq_ch
	main:
		version_ch = Channel.empty()
	
		r1r2_ch= fastq_ch.map{T->T.plus([
			"reference":reference
			])}
	
		bam1_ch = BWA_MEM_01(meta,r1r2_ch)
		version_ch = version_ch.mix(bam1_ch.version)
	
		merge_ch = SAMTOOLS_MERGE_01(meta, reference, bam1_ch.bam.groupTuple())
		version_ch = version_ch.mix(merge_ch.version)

		markdup_ch = SAMBAMBA_MARKDUP_01(meta,merge_ch.bam)
		version_ch = version_ch.mix(markdup_ch.version)

		if( parseBoolean(meta.with_bqsr) ) {
			beds_ch = CONTIGS_IN_BAM(meta,reference,markdup_ch.bam)
			version_ch = version_ch.mix(beds_ch.version)

			brecal_ch = GATK4_BASE_RECALIBRATOR_01(meta, beds_ch.output.splitCsv(header:true,sep:'\t'))
			version_ch = version_ch.mix(brecal_ch.version)
	
			gather_ch = GATK4_GATHER_BQSR_01(meta, brecal_ch.output.map{T->[
				["sample":T[0].sample,"bam":T[0].bam],
				T[1]
				]}.groupTuple())
			version_ch = version_ch.mix(gather_ch.version)

		
			applybqsr_ch = GATK4_APPLY_BQSR_01(meta, gather_ch.output.map{T->
				T[0].plus("table":T[1])
				})
			version_ch = version_ch.mix(applybqsr_ch.version)
		
			cram_ch = SAMTOOLS_BAM_TO_CRAM_01(meta,reference,applybqsr_ch.bam.map{T->[T[0].sample,T[1]]})
			version_ch = version_ch.mix(cram_ch.version)
			}
		else
			{
			cram_ch = SAMTOOLS_BAM_TO_CRAM_01(meta,reference, markdup_ch.bam)
			version_ch = version_ch.mix(cram_ch.version)		
			}

		version_ch = MERGE_VERSION(meta, "Remap", "Remap bam on another reference", version_ch.collect())
	emit:
		version= version_ch.version
		bams = cram_ch.cram
	}

process CONTIGS_IN_BAM {
tag "${sample} ${bam.name}"
executor "local"
input:
	val(meta)
	val(reference)
	tuple val(sample),val(bam)
output:
	path("bam.contigs.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def recalGroup = meta.recalGroup?:"50mb"
"""
hostname 1>&2
${moduleLoad("samtools jvarkit")}
set -o pipefail
mkdir TMP

samtools idxstats "${bam}" |\
	 awk -F '\t' '(\$3!=0 && \$1!="*") {printf("%s\t0\t%s\\n",\$1,\$2);}' |\
	 java -jar \${JVARKIT_DIST}/bedcluster.jar --reference "${reference}" -o TMP --size "${meta.recalGroup?:"50mb"}"
	 
	 echo -e 'sample\tbam\tbed\trecalibration_vcfs\treference' > bam.contigs.tsv
	 find \${PWD}/TMP -type f -name "*.bed" |\
                awk -F '\t' '{printf("${sample}\\t${bam}\\t%s\\t${meta.recalibration_vcfs?:""}\\t${reference}\\n",\$0);}' >> bam.contigs.tsv

###########################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
       	<entry key="name">${task.process}</entry>
        <entry key="description">bam to bed</entry>
        <entry key="recalGroup">${recalGroup}</entry>
</properties>
EOF
"""
}
