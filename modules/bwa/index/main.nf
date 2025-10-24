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
process BWA_INDEX {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/bwa.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(fasta)
output:
	tuple val(meta),path("${fasta.baseName}BWAIndex"),emit:bwa_index
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args = task.ext.args?:""
"""
    mkdir -p TMP
    bwa index ${args} -p TMP/${prefix} ${fasta}
    mv TMP "${fasta.baseName}BWAIndex"

cat << END_VERSIONS > versions.yml
${task.process}:
    bwa : todo
END_VERSIONS
"""

stub:
"""
mkdir -p TMP
touch TMP/${fasta.baseName}.amb
touch TMP/${fasta.baseName}.ann
touch TMP/${fasta.baseName}.bwt
touch TMP/${fasta.baseName}.pac
touch TMP/${fasta.baseName}.sa
mv TMP "${fasta.baseName}BWAIndex"
touch versions.yml
"""
}
