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
array 100
tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(vcf),path(idx)
output:
	tuple val(meta),path("*.bed"),path(vcf),path(idx),optional:true,emit:output
	tuple val(meta),path("*.bed"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:	
	def prefix= task.ext.prefix?:vcf.baseName+".interval"
"""
hostname 1>&2
set -o pipefail
mkdir -p OUT TMP

bcftools index -s "${vcf}" |\\
	awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed
	
if test -s TMP/jeter.bed
then
	mv TMP/jeter.bed "${prefix}.bed"
fi


cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
"""
}
