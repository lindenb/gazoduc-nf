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
include {isBlank    } from '../../utils/functions.nf'
include {verify     } from '../../utils/functions.nf'


process DEPTH_OUTLIERS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.bed.gz"),path("*.tbi"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.outliers"
	def mapq   = task.ext.mapq?:1
    def merge_distance   = task.ext.merge_distance?:0
	def args1  = task.ext.args1?:"-F 3844"
    def awk_expr = task.ext.awk_expr?:""
    verify(!isBlank(awk_expr),"${task.process} task.ext.awk_expr must be defined (e.g \$3>100) .")
"""
hostname 1>&2
mkdir -p TMP

samtools view \\
    --min-MQ "${mapq}" \\
    --uncompressed ${args1} \\
    ${optional_bed?"-L ${optional_bed}":""} \\
    --reference "${fasta}" \\
    "${bam}" |\\
    samtools depth \\
        -q "${mapq}" \\
        ${optional_bed?"-b ${optional_bed}":""}  - |\\
        awk -F '\t' '(${awk_expr} ) {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$2)-1,\$2,\$3);}' |\\
        bedtools merge -d ${merge_distance} -c 4 -o median > TMP/jeter.bed

if test -s TMP/jeter.bed
then
    sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed | bgzip > TMP/jeter.bed.gz
    tabix -f -p bed TMP/jeter.bed.gz
    mv TMP/jeter.bed.gz ${prefix}.bed.gz
    mv TMP/jeter.bed.gz.tbi ${prefix}.bed.gz.tbi
fi


cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
	bedtools: "\$(bedtools --version | awk '{print \$NF}')"
END_VERSIONS
"""

stub:
    def awk_expr = task.ext.awk_expr?:""
    verify(!isBlank(awk_expr),"${task.process} awk_expr must be defined.")
    def prefix = task.ext.prefix?:"${meta.id}.outliers"
"""
touch versions.yml ${prefix}.bed.gz ${prefix}.bed.gz.tbi
"""
}
