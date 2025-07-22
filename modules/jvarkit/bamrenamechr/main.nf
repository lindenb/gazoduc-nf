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

process JVARKIT_BAM_RENAME_CONTIGS{
label "process_short"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
tag "${meta.id?:bam.name}"
input:
	tuple val(meta1),path(dict_source)
	tuple val(meta2),path(optional_bed)
	tuple val(meta ),path(bam),path(bai),path(fasta),path(fai),path(dict)
output:
	tuple val(meta ),path("*.bam"),path("*.bai"),emit:bam
	path("versions.yml"),emit:versions
script:
	def cpus = ((task.cpus?:1) as int)

	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:"-i"
	def prefix  = task.ext.prefix?:bam.baseName+".setdict"
	def has_bed = optional_bed?true:false
"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

	if ${has_bed} 
	then
		jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr \\
					-R ${fasta} -c 1 --convert SKIP  ${optional_bed}  |\\
				cut -f1,2,3 |\\
				sort -T TMP -t '\t' -k1,1 -k2,2n |\\
				bedtools merge > TMP/jeter.bed
	fi


	samtools view ${cpus>1?"--threads ${cpus -1}":""} --uncompressed -T ${fasta} ${args1} ${has_bed?"-M -L TMP/jeter.bed":""} "${bam}" |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bamrenamechr \\
			${args2} \\
			--dict '${dict_source}' \\
			--samoutputformat BAM \\
			-o TMP/jeter.bam

	samtools view -H TMP/jeter.bam > TMP/jeter.header
	
	if head -n1 TMP/jeter.header | grep -F 'SO:coordinate'
	then
		echo "no need to sort" 1>&2
	else
		samtools sort  -m '${task.memory.giga}G' --threads '${task.cpus}' -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
		mv TMP/jeter2.bam TMP/jeter.bam
	fi

	samtools index  -@ ${task.cpus} TMP/jeter.bam

	mv TMP/jeter.bam ${prefix}.bam
	mv TMP/jeter.bam.bai ${prefix}.bam.bai

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
    samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
