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

process GATK4_APPLY_BQSR {
label "process_short"
tag "${row.sample} ${row.bam}"
afterScript 'rm -rf TMP'
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta),path(bam),path(bai),path(table)
output:
    tuple val(row),path("*.bam"),path("*.bai"),emit:bam
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:bam.baseName+".bqsr"
"""
hostname 1>&2
mkdir TMP

# allele specific annotation are not supported in non-gvcf mode
gatk --java-options gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP" ApplyBQSR \
	-R ${fasta} \\
	-I "${row.bam}" \\
	-bqsr "${table}" \\
	-O "TMP/${prefix}.bam"

if test -f "TMP/${prefix}.bai"
then
	mv "TMP/${prefix}.bai" "TMP/${prefix}.bam.bai"
fi

mv "TMP/${prefix}.bam" ./
mv "TMP/${prefix}.bam.bai" ./

cat << EOF > version.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""

}
