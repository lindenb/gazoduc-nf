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
process READ_LENGTH {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
	tuple val(meta),path(bam)
output:
	tuple val(meta),path("${sample}.len.txt"),emit:length
    path("versions.yml"),emit:versions
script:
	def n_reads = task.ext.n_reads?:1000
    def prefix = task.ext.prefix?:"${meta.id}"
"""
samtools view -F 3844 -T "${fasta}" "${bam}" |\\
	cut -f 10 |\\
	head -n "${n_reads}" |\\
	awk -F '\t' 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$1)} END{print (N==0?150:(T/N));}' > ${prefix}.len.txt

touch versions.yml
"""
stub:
     def prefix = task.ext.prefix?:"${meta.id}"
"""
echo 150 >  ${prefix}.len.txt
touch versions.yml
"""
}
