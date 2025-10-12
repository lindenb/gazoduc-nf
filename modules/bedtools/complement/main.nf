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
process BEDTOOLS_COMPLEMENT {
tag "${bed.name}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
     tuple val(meta1),path(fai)
     tuple val(meta ),path(bed)
output:
    tuple val(meta),path("*.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
    def args1 = task.args1?:""
    def prefix = task.ext.prefix?:"${bed.baseName}.complement"
    def awkargs =  task.ext.awkargs?:"(1==1)"
"""
mkdir -p TMP
cut -f1,2 "${fai}" |\\
    sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/sizes.txt

cut -f1,2,3 "${bed}" |\\
sort -T TMP -t '\t' -k1,1 -k2,2n  |\\
	bedtools complement -i - -g TMP/sizes.txt ${args1} |\\
    awk -F '\t' '(int(\$2) < int(\$3) || ${awkargs})' |\\
     sort -T TMP -t '\t' -k1,1 -k2,2n  > TMP/jeter.bed

mv -v TMP/jeter.bed "${prefix}.bed"


cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bedtools --version | awk '(NR==1)  {print \$NF}')
END_VERSIONS
"""

stub:
def prefix = task.ext.prefix?:"${bed.baseName}.complement"
"""
touch "${prefix}.bed"
touch versions.yml
"""
}
