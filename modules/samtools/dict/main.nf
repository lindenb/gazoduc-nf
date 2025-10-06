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
process SAMTOOLS_DICT {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta ),path(fasta)
output:
    tuple val(meta),path("*.dict"),emit:dict
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:"-A"
    def args2 = task.ext.args2?:(meta.uscsc_name?"-a ${meta.uscsc_name}":"")
"""

samtools dict ${args1} ${args2} -o "${fasta.baseName}.dict" ${fasta}

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
"""
touch "${fasta.baseName}.dict"
touch versions.yml
"""
}