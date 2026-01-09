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

include {verify } from '../../utils/functions.nf'

process BCFTOOLS_VIEW {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(bed)
    tuple val(meta2),path(samples)
	tuple val(meta ),path(vcf),path(opt_idx)
output:
	tuple val(meta),path("*.vcf.gz"),emit:vcf
    tuple val(meta),path("*.vcf.gz.tbi"),optional:true,emit:tbi
	path("versions.yml"),emit:versions
script:
    def negate_samples = task.ext.negate_samples?:"" // or use '^'
    verify(negate_samples=="" || negate_samples=="^", "${task.process} negate_samples expected empty or '^'")
    def negate_bed = task.ext.negate_bed?:""// or use '^'
    verify(negate_bed=="" || negate_bed=="^", "${task.process} negate_bed expected empty or '^'")
    def bed_arg  = task.ext.bed_arg?:"regions"
	def prefix = task.ext.prefix?:"${meta.id}.view"
    def args1 = task.ext.args1?:""
    def args2 = task.ext.args2?:""
    def args3 = task.ext.args3?:""
"""
mkdir -p TMP
set -x

bcftools view \\
    --threads ${task.cpus} \\
    ${args1} \\
    ${args2} \\
    ${args3} \\
    ${samples?"--samples-file \"${negate_samples}${samples}\"":""} \\
    ${bed?"--${bed_arg}-file \"${negate_bed}${bed}\"":""} \\
    -O z \\
    -o TMP/jeter.vcf.gz \\
    "${vcf}"

mv  TMP/jeter.vcf.gz ${prefix}.vcf.gz

if test -f TMP/jeter.vcf.gz.tbi
then
    mv  TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi
fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
"""
touch versions.yml "${meta.id}.vcf.gz"
"""
}
