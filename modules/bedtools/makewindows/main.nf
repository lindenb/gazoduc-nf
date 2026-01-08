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
process BEDTOOLS_MAKEWINDOWS {
tag "${bed.name}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(bed)
output:
    tuple val(meta),path("*.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
    def option1= (bed.name.endsWith(".fai")?"-g":"-b")
    def args = task.ext.args?:""
    if((args as String).trim().isEmpty()) throw new IllegalArgumentException("args empty for ${task.process}")
    def prefix = task.ext.prefix?:"${meta.id?:bed.baseName}.makewindows"
"""
mkdir -p TMP

bedtools makewindows ${args} ${option1} "${bed}" |\\
    sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed

mv TMP/jeter.bed "${prefix}.bed"


cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bedtools --version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
"""
cp "${bed}" "${meta.id?:bed.baseName}.makeWindows.bed"
touch versions.yml
"""
}
