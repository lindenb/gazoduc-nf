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


process BEDTOOLS_INTERSECT {
label "process_single"
tag "${meta.id?:""} ${file1.name} ${file2.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(optional_fai)
	tuple val(meta),path(file1),path(file2)    
output:
	tuple val(meta),path("*.bed"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:	
    def sizes = optional_fai?"-g ${optional_fai}":''
    def prefix = task.ext.prefix?:"\${MD5}.intersect"
    def args1 = task.ext.args1?:""
"""
hostname 1>&2
mkdir -p TMP

bedtools intersect \\
    ${sizes} \\
    ${args1} \\
    -a ${file1} \\
    -b ${file2} > "TMP/jeter.bed"


if test -s TMP/jeter.bed
then

sort -T TMP -k1,1 -k2,2n -t '\t' TMP/jeter.bed > TMP/jeter2.bed

MD5=`cat TMP/jeter2.bed | md5sum | cut -d ' ' -f1`


mv TMP/jeter2.bed "${prefix}.bed"

fi

cat << EOF > versions.yml
${task.process}:
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
EOF
"""

stub:
def prefix="intersection"
"""
cp ${file1} "${prefix}.bed"
touch versions.yml
"""
}
