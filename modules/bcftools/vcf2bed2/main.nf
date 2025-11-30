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


process VCF_TO_BED {
label "process_single"

tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(vcf)
output:
	tuple val(meta),path("*.bed"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:	
    def args1 = task.ext.args1?:""
    def args2 = task.ext.args1?:""
	def prefix= task.ext.prefix?:"${meta.id}.vcf2bed"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP


bcftools query ${args1} -f '%CHROM\t%POS0\t%END\\n' "${vcf}" |\\
	LC_ALL=C sort --buffer-size=${task.memory.mega}M -T TMP -t '\t' -k1,1 -k2,2n |\\
    bedtools merge > TMP/jeter.bed


if test -s TMP/jeter.bed
then
		mv TMP/jeter.bed "${prefix}.bed"
fi


cat << EOF > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
EOF
"""
stub:
def prefix= task.ext.prefix?:"${meta.id}.vcf2bed"
"""
touch versions.yml ${prefix}.${meta.id}.bed
"""
}
