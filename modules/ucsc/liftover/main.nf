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

process LIFTOVER {
tag "${bed.name} ${chain.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(chain)
	tuple val(meta ),path(bed)
output:
	tuple val(meta),path("liftover.*"),emit:bed
	tuple val(meta),path("unmapped.*"),optional:true,emit:unmapped
	path("versions.yml"), emit:versions
script:
	def args1 = task.ext.args1?:""
	def suffix = "${meta.id}.lift${meta1.id}"
"""
hostname 1>&2
mkdir TMP
liftOver ${args1} "${bed}" "${chain}" TMP/liftover.bed TMP/unmapped.bed


sort -T TMP -t '\t' -k1,1 -k2,2n TMP/liftover.bed > TMP/jeter.bed
mv TMP/jeter.bed liftover.${suffix}.bed

sort -T TMP -t '\t' -k1,1 -k2,2n TMP/unmapped.bed > TMP/jeter.bed
mv TMP/jeter.bed unmapped.${suffix}.bed



cat << EOF > versions.yml
"${task.process}":
	ucsc
EOF
"""

stub:
	def suffix = "${meta.id}.lift${meta1.id}"
"""
touch unmapped.${suffix}.bed remapped.${suffix}.bed  versions.yml
"""
}
