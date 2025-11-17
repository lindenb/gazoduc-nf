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

process BAM_TO_CRAM {
tag "${meta.id} ${bam.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta ),path(bam)
output:
	tuple val(meta),path("*.cram"), path("*.crai"),emit:cram
	tuple val(meta),path("*.md5"),emit:md5
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def level = task.ext.compression_level?:9
	def prefix = task.ext.prefix?:meta.id
"""
hostname 1>&2
mkdir -p TMP
samtools view \\
	${args1} \\
	-@ ${task.cpus} \\
	--write-index \\
	-O "CRAM,level=${level}" \\
	-o TMP/jeter.cram \\
	-T "${fasta}" \\
	"${bam}"


mv TMP/jeter.cram "${prefix}.cram"
mv TMP/jeter.cram.crai "${prefix}.cram.crai"
md5sum "${prefix}.cram" > "${prefix}.md5"

cat << EOF > versions.yml
"${task.process}":
    samtools:\$(samtools  --version | head -n 1| cut -d ' ' -f2)
EOF
"""

stub:
"""
touch ${meta.id}.md5
touch ${meta.id}.cram
touch ${meta.id}.cram.crai
touch versions.yml
"""
}

