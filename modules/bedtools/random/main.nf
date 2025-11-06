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


process BEDTOOLS_RANDOM {
label "process_single"
tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(fai)
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:	
    def prefix = task.ext.prefix?:"${meta.id}"
    def args1 = task.ext.args1?:""
"""
mkdir -p TMP

bedtools random \\
    ${args1} \\
    -g ${fai} |\\
    sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed

mv TMP/jeter.bed "${prefix}.bed"

cat << EOF > versions.yml
${task.process}:
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch "${prefix}.bed" versions.yml
"""
}
