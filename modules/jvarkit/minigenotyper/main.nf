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
include {verify           } from '../../../modules/utils/functions.nf'

process MINIGENOTYPER {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(vcf)
    tuple val(meta),path("BAMS/*")
output:
    tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.${meta4.id}"
    def jvm= task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
	def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"
    def args1 = task.ext.args1?:"--skip-illegal-variant --maxRecordsInRam 100000 "
    def args2 = task.ext.args2?:""
"""

mkdir -p TMP
find BAMS/ -name "*am" | sort -V > TMP/bams.list
test -s TMP/bams.list

 ${jvarkit} minigenotyper \\
    -V ${vcf} \\
    -R ${fasta} \\
    ${args1} \\
    TMP/bams.list  > TMP/jeter.vcf

bcftools sort ${args2} --max-mem ${task.memory.mega}M -T TMP/sort -O z -o TMP/jeter2.vcf.gz TMP/jeter.vcf
bcftools index --threads ${task.cpus} -f -t  TMP/jeter2.vcf.gz


mv TMP/jeter2.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter2.vcf.gz.tbi "${prefix}.vcf.gz.tbi"

cat << END_VERSIONS > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
END_VERSIONS
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.${meta4.id}"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}
