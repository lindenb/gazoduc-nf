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


process SAMTOOLS_STATS {
label "process_single"
tag "${meta.id?:bam.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.stats.tsv"),emit:stats
	path("versions.yml"),emit:versions
script:
	def args = task.ext.args?:""
	def prefix = meta.id?:bam.baseName
	def has_bed=optional_bed?true:false
"""
hostname 1>&2
mkdir -p TMP
if ${has_bed}
then
	# ONE based intervals ...
	${has_bed && optional_bed.name.endsWith(".gz")?"gunzip -c":"cat"} ${optional_bed} |\\
		awk -F '\t' '{printf("%s\t%s\t%s\\n",\$1,int(\$2)+1,\$3);}' > TMP/jeter.intervals
fi

samtools stats \\
	${args} \\
	${has_bed?"--target-regions TMP/jeter.intervals":""} \\
	--ref-seq ${fasta} \\
	--reference ${fasta} \\
	--threads ${task.cpus} "${bam}" > TMP/jeter.tsv

mv TMP/jeter.tsv ${prefix}.stats.tsv


cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
"""
touch versions.yml ${meta.id}.stats.tsv
"""
}
