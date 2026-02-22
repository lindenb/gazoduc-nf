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

process VARIANT_RECALIBRATOR {
tag "${meta.id} ${meta.type}"
label "process_short"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(other_arguments)
	tuple val(meta),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path("*.tranches.txt"),emit:vcf
	tuple val(meta),path("*.plot.R"),emit:R
	path("versions.yml"),emit:versions
script:
    def args1= task.ext.args1?:""
	def type = task.ext.type?:(meta.type?:"")
	verify(!isBlank(type),"undefined meta.type in ${task.process}")
    verify(type.matches("(snps|indels)"),"${type} mus be snps or indels in ${task.process}")
	def prefix = task.ext.prefix?:"${meta.id}.${type}.recal"
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
"""
hostname 1>&2
set -x
mkdir -p TMP

gatk  --java-options "${jvm}" VariantRecalibrator  \\
	-R "${fasta}" \\
	-V "${vcf}" \\
    ${args1} \\
	${other_arguments?"--arguments_file ${other_arguments}":""} \\
	-mode ${type.equals("snps")?"SNP":"INDEL"} \\
	-O "${prefix}.vcf.gz" \\
	--tranches-file "${prefix}.tranches.txt" \\
    --dont-run-rscript \\
	--rscript-file  "${prefix}.plot.R"

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""
stub:
def prefix = task.ext.prefix?:vcf.name.md5()+".snp.recal" 
"""
touch versions.yml  minikit.jar ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi ${prefix}.plot.R ${prefix}.tranches.txt
"""
}