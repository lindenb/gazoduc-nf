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
include { removeCommonSuffixes } from '../../utils/functions.nf'

process BEDTOOLS_MULTIINTER {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(opt_fai)
    tuple val(meta ),path(beds)
output:
	tuple val(meta),path("*.bed.gz"),path("*.tbi"),emit:bed
	path("versions.yml"),emit:versions
script:
	def L= (beds instanceof List?beds:[beds])
	def L2 = L.sort{f1,f2->f1.name<=>f2.name}
	def args = task.ext.args?:""
	def prefix = task.ext.prefix?:"${meta.id}.multiinter"
	def with_name = (task.ext.with_name?:false)
	def names0 = L2.collect{f->removeCommonSuffixes(f.name)}.join(" ")
	def names = (task.ext.names?:names0)
"""
mkdir -p TMP


bedtools multiinter \\
	${args} \\
	${opt_fai?"-g \"${opt_fai}\" -empty ":""} \\
	${with_name?"-names ${names}":""} \\
	-i ${L2.collect{f->f.name}.join(" ")} |\\
		bgzip > TMP/jeter.bed.gz

tabix -f -p bed TMP/jeter.bed.gz

mv TMP/jeter.bed.gz "${prefix}.bed.gz"
mv TMP/jeter.bed.gz.tbi "${prefix}.bed.gz.tbi"

cat << EOF > versions.yml
${task.process}:
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
EOF
"""


stub:
"""
touch versions.yml
touch ${meta.id}.${meta.treshold}.bed.gz
touch ${meta.id}.${meta.treshold}.bed.gz.tbi
"""
}
