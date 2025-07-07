// cf. workspace/pyflow.data/state/pyflow_tasks_info.txt
workflow {
	samplesheet_ch=Channel.fromPath(params.samplesheet).
		splitCsv(header:true,sep:'\t')


	SINGLE_SAMPLE(
		file(params.fasta),
		file(params.fai),
		file(params.bed),
		samplesheet_ch.map{[it.sample,it.bam,it.bai]}
		)
	}


workflow SINGLE_SAMPLE {
	take:
		fasta
		fai
		bed
		sample_bam_ch
	main:
		aln_ch = ALIGNMENT_STATS(fasta,fai,sample_bam_ch)
		merged_aln_ch = MERGE_ALIGNMENT_STATS(aln_ch.output.groupTuple())
		summarize_ch = SUMMARIZE_STATS(merged_aln_ch.output)
		faidx_ch = FAIDX_INFO(fai)
		contigs_ch  = faidx_ch.output.splitCsv(header:false,sep:'\t')

		chrom_depth = ESTIMATE_CHROM_DEPTH(
			fasta,
			fai,
			sample_bam_ch.join(summarize_ch.output).combine(contigs_ch)
			)
		CAT_CHROM_DEPTH(chrom_depth.output.groupTuple())

	}


process ALIGNMENT_STATS {
tag "${sample}"
label "process_single"
conda "${moduleDir}/../../conda/manta.yml"
input:
	path(fasta)	
	path(fai)	
	tuple val(sample),path(bam),path(bai)
output:
	tuple val(sample),path("${sample}.alignmentStats.xml"),emit:output
script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../libexec/GetAlignmentStats \\
	--ref "${fasta}" \\
	--output-file "${sample}.alignmentStats.xml" \\
	--align-file "${bam}"
"""
}


process MERGE_ALIGNMENT_STATS {
tag "${sample}"
label "process_single"
conda "${moduleDir}/../../conda/manta.yml"
input:
	tuple val(sample),path("STATS/*")
output:
	tuple val(sample),path("${sample}.merged.alignmentStats.xml"),emit:output
script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../libexec/MergeAlignmentStats \\
	`find ./STATS -name "*.xml" -printf "--align-stats-file  %p "` \\
	--output-file "${sample}.merged.alignmentStats.xml"
"""
}


process SUMMARIZE_STATS {
tag "${sample}"
label "process_single"
conda "${moduleDir}/../../conda/manta.yml"
input:
        tuple val(sample),path(stats)
output:
        tuple val(sample),path("${sample}.alignmentStatsSummary.txt"),emit:output
script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../libexec/SummarizeAlignmentStats \\
	--align-stats "${stats}" \\
	--output-file "${sample}.alignmentStatsSummary.txt"
"""
}

process FAIDX_INFO {
tag "${fai}"
label "process_single"
input:
        path(fai)
output:
        path("contigs.tsv"),emit:output
script:
"""
awk -F '\t' '(\$1 ~ /^(chr)?[0-9XY]*\$/ ) {printf("%s\t%d\\n",\$1,NR);}' '${fai}' > contigs.tsv
"""
}

process ESTIMATE_CHROM_DEPTH {
tag "${sample} ${contig}"
label "process_single"
conda "${moduleDir}/../../conda/manta.yml"
input:
        path(fasta)
        path(fai)
        tuple val(sample),path(bam),path(bai),path(stats),val(contig),val(tid)
output:
        tuple val(sample),path("${sample}.chromDepth.txt.000.txt_${tid}_${contig}.txt"),emit:output
script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../libexec/GetChromDepth \\
	--ref "${fasta}" \\
	--align-file "${bam}" \\
	--output "${sample}.chromDepth.txt.000.txt_${tid}_${contig}.txt" \\
	--chrom "${contig}"
"""
}

process CAT_CHROM_DEPTH {
tag "${sample}"
label "process_single"
conda "${moduleDir}/../../conda/manta.yml"
input:
        tuple val(sample),path("CONTIG/*")
output:
        tuple val(sample),path("${sample}.chromDepth.merge.txt"),emit:output
script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../libexec/cat.py \\
        --output "${sample}.chromDepth.merge.txt" \\
	CONTIG/*.txt
        
"""
}


makeGraph {

script:
"""
CONFIG_MANTA_EXE=`which configManta.py`
MANTA_BIN_DIR=`dirname \${CONFIG_MANTA_EXE}`

\${MANTA_BIN_DIR}/../EstimateSVLoci \\
	--output-file svLocusGraph.chromId_000_chr1_0000.bin \\
	--align-stats ${start} \\
	--region ${interval}  \
	--min-candidate-sv-size 8 \\
	--min-edge-observations 3 \\
	--ref "${fasta}" \\
	--align-file "${bam}" \\
	--chrom-depth "${depth"}
"""
}
