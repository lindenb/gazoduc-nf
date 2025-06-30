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

process BCFTOOLS_BCSQ {
tag "${meta.id}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
   	tuple val(meta3),path(gff3),path(gff3_idx) 
    tuple val(meta ),path(vcf),path(tbi)

output:
    tuple val(meta ),path("*.bcf") ,path("*.csi"),emit:vcf
script:
    def prefix=task.ext.prefix?:vcf.baseName+".bcsq"
"""
hostname 1>&2
mkdir -p TMP
              
bcftools csq \\
    --threads ${task.cpus} \\
    -O b \\
    --force \\
    --local-csq \\
    --ncsq 10000 \\
    --fasta-ref "${fasta}" \\
    --gff-annot "${gff3}" \\
    -o TMP/${prefix}.bcf \\
    "${vcf}"

bcftools index \\
    --threads ${task.cpus} \\
    --force TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./
"""
}
