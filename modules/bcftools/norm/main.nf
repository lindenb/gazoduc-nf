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
process BCFTOOLS_NORM {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.bcf"),path("*.bcf.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:"--remove-duplicates --multiallelics -both --old-rec-tag MULTIALLELIC"
	def args3 = task.ext.args3?:"-i 'ALT!=\"*\"'"
	def prefix = task.ext.prefix?:vcf.baseName+".norm"
	def argsbed = optional_bed?"--regions-file \"${optional_bed}\"":""
	def set_id = task.ext.set_id?:"" //eg %VKX
"""
mkdir -p TMP
set -o pipefail

bcftools view ${argsbed} ${args1} -O u '${vcf}' |\\
	bcftools norm ${args2} --fasta-ref '${fasta}'  -O u |\\
	${set_id.isEmpty()?"":"bcftools annotate --set-id '${set_id}' -O u |"}\\
	bcftools view ${args3} --write-index  -O b -o TMP/${prefix}.bcf

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
