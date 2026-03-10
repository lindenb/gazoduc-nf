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

include { verify             } from '../../../modules/utils/functions.nf'
include { isBlank            } from '../../../modules/utils/functions.nf'


process APPLY_VQSR {
tag "${meta.id} ${meta4.type?:""}"
label "process_single"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
    tuple val(meta4),path(recal_vcf),path(recal_vcf_idx),path(tranches)
	tuple val(meta), path(vcf),path(vcfidx)
output:
	tuple val(meta), path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def type = task.ext.type?:(meta4.type?:"")
	verify(!isBlank(type),"undefined meta4.type in ${task.process}")
    verify(type.matches("(snps|indels)"),"${type} mus be snps or indels in ${task.process}")
	def prefix = task.ext.prefix?:"${meta.id}.${type}.vqsr"
    def level = task.ext.level?:99.0
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
"""
hostname 1>&2
mkdir -p TMP


gatk  --java-options "${jvm}" ApplyVQSR  \\
        -R "${fasta}" \\
        -V ${vcf} \\
        ${args1} \\
        -mode ${type.equals("snps")?"SNP":"INDEL"} \\
        --truth-sensitivity-filter-level ${level} \\
        --recal-file "${recal_vcf}" \\
        --tranches-file "${tranches}" \\
        -O "TMP/jeter2.vcf.gz"

mv TMP/jeter2.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi ${prefix}.vcf.gz.tbi

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( (gatk --java-options "${jvm}" --version 2> /dev/null  | paste -s -d ' ' ) || true)"
EOF
"""

stub:
    def type = task.ext.type?:(meta4.type?:"")
	verify(!isBlank(type),"undefined meta4.type in ${task.process}")
    verify(type.matches("(snps|indels)"),"${type} mus be snps or indels in ${task.process}")
	def prefix = task.ext.prefix?:"${meta.id}.${type}.vqsr"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}