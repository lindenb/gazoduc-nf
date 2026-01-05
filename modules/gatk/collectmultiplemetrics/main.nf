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
process COLLECT_MULTIPLE_METRICS {
tag "${meta.id?:""}"
label "process_short" // longer than 3H for a wgs
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(optional_refflat)
    tuple val(meta5),path(optional_intervals)
    tuple val(meta),path(bam),path(bai)
output:
	tuple val(meta),path("*.alignment_summary_metrics"),optional:true,emit:alignment_summary_metrics
	tuple val(meta),path("*.base_distribution_by_cycle_metrics"),optional:true,emit:base_distribution_by_cycle_metrics
	tuple val(meta),path("*.base_distribution_by_cycle.pdf"),optional:true,emit:base_distribution_by_cycle_pdf
	tuple val(meta),path("*.insert_size_histogram.pdf"),optional:true,emit:insert_size_histogram_pdf
	tuple val(meta),path("*.insert_size_metrics"),optional:true,emit:insert_size_metrics
	tuple val(meta),path("*.quality_by_cycle_metrics"),optional:true,emit:quality_by_cycle_metrics
	tuple val(meta),path("*.quality_by_cycle.pdf"),optional:true,emit:quality_by_cycle_pdf
	tuple val(meta),path("*.quality_distribution_metrics"),optional:true,emit:quality_distribution_metrics
	tuple val(meta),path("*.quality_distribution.pdf"),optional:true,emit:quality_distribution_pdf
	tuple val(meta),path("*.read_length_histogram.pdf"),optional:true,emit:read_length_histogram_pdf
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
    def args1 = task.ext.args1?:"--PROGRAM CollectAlignmentSummaryMetrics --PROGRAM CollectBaseDistributionByCycle --PROGRAM CollectInsertSizeMetrics --PROGRAM MeanQualityByCycle --PROGRAM QualityScoreDistribution"
    def args2 = task.ext.args2?:""
"""
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP" CollectMultipleMetrics \\
    --REFERENCE_SEQUENCE "${fasta}" \\
    --INPUT "${bam}" \\
    --OUTPUT "${prefix}" \\
    ${optional_refflat?"--REF_FLAT \"${optional_refflat}\"":""} \\
    ${optional_intervals?"--INTERVALS \"${optional_intervals}\"":""} \\
    ${args1} \\
    ${args2} 1>&2

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""


stub:
"""
touch versions.yml
touch ${meta.id}.alignment_summary_metrics
touch ${meta.id}.base_distribution_by_cycle_metrics
touch ${meta.id}.insert_size_metrics
touch ${meta.id}.quality_by_cycle_metrics
touch ${meta.id}.quality_distribution_metrics
"""
}
