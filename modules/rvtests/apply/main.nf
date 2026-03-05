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
process RVTESTS_APPLY {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/rvtests.yml"
input:
    tuple val(meta1), path(pedigree)
    tuple val(meta2), path(vcf),path(idx)
    tuple val(meta ), path(setfile)
output:
    tuple val(meta),path("*.assoc", arity: '1..*'),emit:assoc
    path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def prefix= task.ext.prefix?:"${meta.id}.rvtests"
"""

rvtest  --noweb \
    --inVcf "${vcf}" \
    --setFile "${setfile}" \
    --pheno "${pedigree}" \
    --out "${prefix}." \
	${args1} 2> ${prefix}.rvtest.log

cat << EOF > versions.yml
${task.process}:
    rvtest: todo
EOF
"""
stub:
    def prefix= task.ext.prefix?:"${meta.id}.rvtests"
"""
touch versions.yml "${prefix}.assoc"
"""
}