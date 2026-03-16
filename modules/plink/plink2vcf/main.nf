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
 
process PLINK_RECODE_VCF {
label "process_single"
afterScript "rm -rf TMP"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta ),path(bim),path(bed),path(fam)
output:
    tuple val(meta),path("*.vcf.gz"),emit:vcf
    tuple val(meta),path("*.log"),emit:log
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.arg1?:""
    def args2 = task.ext.arg2?:" vcf-iid  "
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

plink  \\
	${args1} \\
    --real-ref-alleles \\
    --bfile ${bim.baseName} \\
	--recode ${args2} bgz \\
    --out TMP/jeter \\
	--threads ${task.cpus} > ${prefix}.log

mv TMP/jeter.vcf.gz ${prefix}.vcf.gz

cat << EOF > versions.yml
${task.process}:
    plink: \$(plink --version)
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.vcf.gz 
"""
}
