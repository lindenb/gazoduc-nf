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


process SAMTOOLS_DEPTH{
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.depth.tsv"),emit:depth
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
	def mapq   = task.ext.mapq?:1
	def args1  = task.ext.args1?:"-F 3844 -M"
	def args2  = task.ext.args2?:"-a"
"""
hostname 1>&2

if ${optional_bed?true:false}
then

	samtools view ${args1} --min-MQ "${mapq}" --uncompressed -O BAM  -L ${optional_bed} --reference "${fasta}" "${bam}" |\\
		samtools depth ${args2} -q "${mapq}" -b "${optional_bed}" -  > "${prefix}.depth.tsv"

else

	samtools depth ${args2} -q "${mapq}" "${bam}"  > "${prefix}.depth.tsv"

fi


cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
