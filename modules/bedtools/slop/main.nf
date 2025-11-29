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
include {verify        } from '../../utils/functions.nf'

process BEDTOOLS_SLOP {
label "process_single"
tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fai)    
	tuple val(meta ),path(bed)    
output:
	tuple val(meta),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:	
    def args = task.ext.args?:""
    verify(task.ext.slop!=null,"${task.process} ext.slop cannot be blank e.g:  10")
    def slop = task.ext.slop?:0
    def prefix = task.ext.prefix?:"${meta.id}.slop${slop}"

"""
hostname 1>&2
mkdir -p TMP

cut -f1,2 "${fai}" |\\
    sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.genome

bedtools slop -i "${bed}" -g TMP/jeter.genome  ${args} -b ${slop} > TMP/jeter.bed


mv TMP/jeter.bed "${prefix}.bed"


cat << EOF > versions.yml
${task.process}:
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
EOF
"""

stub:
    verify(task.ext.slop!=null,"${task.process} ext.slop cannot be blank e.g:  10")
    def prefix="slop"
"""
touch "${prefix}.bed"  versions.yml
"""
}
