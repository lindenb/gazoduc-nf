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

process SAM_TO_FASTQ {
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
label "process_short"
afterScript 'rm -rf TMP'
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),
			path("*.R1.fq.gz"),
			path("*.R2.fq.gz"),
            path("*.R0.fq.gz"),//unpaired
            optional:true,emit:paired_end
    tuple val(meta),
			path("*.SE.fq.gz"),
            optional:true,emit:single_end
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.bam2fq"
	def args1 =  task.ext.args1?:""
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

set +o pipefail
samtools view \\
        --threads ${task.cpus} \\
        ${fasta?"--reference \"${fasta}\"":""} \\
        --require-flags 1 \\
        --uncompressed \\
        ${bam} |\\
    samtools head -h 0 -n 1 > TMP/paired.flag

set -o pipefail

if test -s TMP/paired.flag

    gatk --java-options "${jvm}" SamToFastq \\
        ${args1} \\
        --INPUT "${bam}" \\
        --FASTQ "TMP/${prefix}.R1.fq.gz" \\
        ${fasta?"REFERENCE_SEQUENCE \"${fasta}\"":""} \\
        --SECOND_END_FASTQ "TMP/${prefix}.R2.fq.gz" \\
        --UNPAIRED_FASTQ "TMP/${prefix}.R0.fq.gz\" \\
        --VALIDATION_STRINGENCY LENIENT

else

    gatk --java-options "${jvm}" SamToFastq \\
        ${args1} \\
        --INPUT "${bam}" \\
        --FASTQ "TMP/${prefix}.SE.fq.gz" \\
        ${fasta?"REFERENCE_SEQUENCE \"${fasta}\"":""} \\
        --VALIDATION_STRINGENCY LENIENT

fi



mv TMP/*.fq.gz ./ 

 
cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""

stub:
 def prefix = task.ext.prefix?:"${meta.id}.bam2fq"
 def paired = ((task.ext.paired?:true) as boolean)
"""
if ${paired}
then
    touch "${prefix}.R1.fq.gz"
    touch "${prefix}.R2.fq.gz"
    touch "${prefix}.R0.fq.gz"
else
    touch "$prefix}.fq.gz"
fi

touch versions.yml
"""
}
