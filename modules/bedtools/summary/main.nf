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


process BEDTOOLS_SUMMARY {
label "process_single"
tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fai)
	tuple val(meta),path(bed)    
output:
	tuple val(meta),path("*.tsv"),emit:summary
	path("versions.yml"),emit:versions
script:	
    def prefix = task.ext.prefix?:bed.baseName+".summary"
    def awk_expr = task.ext.awk_expr:?""
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

awk '${awk_expr}{printf("%s\t%s\\n",\$1,\$2);}' > TMP/genome.tsv

bedtools summary \\
    -i ${bed} \\
    -g TMP/genome.tsv > "${prefix}.tsv"


cat << EOF > versions.yml
${task.process}:
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
EOF
"""
}
