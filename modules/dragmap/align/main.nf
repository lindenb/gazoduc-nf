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
process DRAGMAP_ALIGN {
label "process_short"
conda "${moduleDir}/../../../conda/dragmap.yml"
tag "${meta.id} ${R1.name} ${R2?R2.name:""}"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(index_dir)
	tuple val(meta4),path(optional_bed)
	tuple val(meta ),path(R1),path(R2)
	
output:
	tuple val(meta),path("*.bam"), path("*.bai"),emit:bam
	path("versions.yml"),emit:versions
script:
	def ID = task.ext.ID?:"${meta.id}"
	def CN = task.ext.CN?:(meta.CN?:"Nantes")
	def PL = task.ext.PL?:(meta.PL?:"ILLUMINA")
	def LB = task.ext.LB?:(meta.LB?:ID)
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args1?:""
	def cpus1 = ((task.cpus?:2) as int)
	def cpus2 = cpus1 -1
	def prefix = task.ext.prefix?:"${ID}${optional_bed?".${optional_bed.name}":""}"
"""
hostname 1>&2
mkdir -p TMP
set -x


dragen-os \\
    --RGID ${ID} \\
    --RGSM ${ID} \\
    --num-thread ${cpus1} \\
    -r ${index_dir} \\
    -1 "${R1}" \\
    ${R2?"-2 ${R2}":""} |\\
	samtools view -O BAM ${args2} ${optional_bed?"-L ${optional_bed}":""} -o TMP/jeter.bam


samtools sort --threads ${task.cpus} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
mv TMP/jeter2.bam TMP/jeter.bam

samtools index  -@ ${task.cpus} TMP/jeter.bam

mv "TMP/jeter.bam" "${prefix}.sorted.bam"
mv "TMP/jeter.bam.bai" "${prefix}.sorted.bam.bai"

##################
cat << END_VERSIONS > versions.yml
${task.process}:
    dragmap : \$(dragen-os -V)
END_VERSIONS
"""

stub:
"""
touch "${meta.id}.sorted.bam"
touch "${meta.id}.sorted.bam.bai"
touch versions.yml
"""
}

