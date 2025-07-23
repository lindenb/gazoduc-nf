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

process BCFTOOLS_STATS {
label "process_single"
tag "${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta4),path(optional_gtf)
	tuple val(meta5),path(optional_samples)
	tuple val(meta ),path("VCFS/*")
output:
    tuple val(meta),path("*.txt"),emit:stats
	path("versions.yml"),emit:versions
script:
	def with_gtf = optional_gtf?true:false
	def with_samples = optional_samples?true:false
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:""
	def prefix = task.ext.prefix?:meta.id+".stats"
	
"""
	hostname 1>&2
	mkdir TMP

	find VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list

	if ${with_gtf}
	then
		gunzip -c ${optional_gtf} |\\
			awk -F '\t' '(\$3=="exon") {printf("%s\t%s\t%s\\n",\$1,\$4,\$5);}' |\\
			sort -T TMP -k1,1 -k2,2n |\\
			bedtools merge > TMP/exons.bed
	fi


	bcftools concat ${args1} -O u -a \\
		--threads ${task.cpus} \\
		--remove-duplicates \\
		${optional_bed?"--regions-file ${optional_bed}":""} \\
		--file-list  TMP/jeter.list |\\
	bcftools stats \\
		${args2} \\
		--threads ${task.cpus} \\
		${with_samples?"--samples-file ${optional_samples}":""} \\
		${with_gtf?"--exons TMP/exons.bed":""} \\
		-F ${fasta} > 'TMP/${prefix}.txt'

mv TMP/${prefix}.txt ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
