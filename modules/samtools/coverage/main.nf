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


process SAMTOOLS_COVERAGE {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta ),path("BAMS/*")
output:
	tuple val(meta),path("*.txt"),emit:output
	path("versions.yml"),emit:versions
script:
	def interval = task.ext.interval?:"${meta.interval?:""}"
	def prefix = task.ext.prefix?:"${meta.id}"
	def mapq   = task.ext.mapq?:0
	def args1  = task.ext.args1?:""
	def width = task.ext.width?:"100"
"""
hostname 1>&2
mkdir -p BAMS TMP
find BAMS/ \\( -name "*.bam" -o -name "*.cram" \\) | sort -T . -V | while read F
do

echo "\${F}" |\\
	samtools samples |\\
	awk -F '\t' '{printf("# %s ${interval} %s\\n",\$1,\$2);}' >> TMP/jeter.txt

samtools coverage \\
	${args1} \\
	--reference "${fasta}" \\
	--min-MQ '${mapq}' \\
	--n-bins '${width}' \\
	${!interval.isEmpty()?"--region \"${interval}\"":""} "\${F}" >> TMP/jeter.txt

echo >> TMP/jeter.txt

done

mv TMP/jeter.txt "${prefix}.txt"

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch "${prefix}.txt" versions.yml
"""
}
