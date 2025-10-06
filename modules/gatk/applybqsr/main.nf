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

process APPLY_BQSR {
label "process_short"
tag "${meta.id} ${bam.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript 'rm -rf TMP'
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta),path(bam),path(bai),path(table)
output:
    tuple val(meta),path("*.bam"),path("*.bai"),emit:bam
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
	def prefix = task.ext.prefix?:bam.baseName+".bqsr"
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

gatk --java-options "${jvm}" ApplyBQSR \
	-R "${fasta}" \\
	-I "${bam}" \\
	-bqsr "${table}" \\
	-O "TMP/${prefix}.bam"

samtools index --threads ${task.cpus} TMP/${prefix}.bam

mv -v "TMP/${prefix}.bam" ./
mv -v "TMP/${prefix}.bam.bai" "${prefix}.bam.bai"

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""



stub:
"""
touch "${meta.prefix}.bqsr.bam"
touch "${meta.prefix}.bqsr.bam.bai"
touch versions.yml
"""
}
