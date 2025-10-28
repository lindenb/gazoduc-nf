/*

Copyright (c) 2024 Pierre Lindenbaum

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


process BAM_TO_FASTQ {
tag "${meta.id?:bam.name}"
label "process_medium"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(bam),path(bai),path(fasta),path(fai),path(optional_bed)
output:
	tuple val(meta),
			path("*.paired.R1.fq.gz"),
			path("*.paired.R2.fq.gz"),
			path("*.singleton.fq.gz"),
			path("*.other.fq.gz"),emit:fastq
	path("versions.yml"),emit:versions
script:
	def has_bed = optional_bed?true:false
	def reference_arg = (fasta?"--reference \"${fasta}\"":"")
	def fetch_pairs = task.attempt==1?" --fetch-pairs ":""
	def prefix = task.ext.prefix?:(meta.id?:bam.baseName)
"""
hostname 1>&2
set -o pipefail

mkdir -p TMP


# extract reads
if ${has_bed}
then
	samtools view  \\
		--fast \\
		--threads ${task.cpus} \\
		${fetch_pairs} \\
		${reference_arg} \\
		-M \\
		-O BAM \\
		-o TMP/jeter.bam \\
		--regions-file "${optional_bed}" \\
		"${bam}"
fi


# collate
samtools collate \\
	-f \\
	--threads ${(task.cpus as int) -1} \\
	-O \\
	-u -\\
	-no-PG \\
	${reference_arg} \\
	"${has_bed ? "TMP/jeter.bam":"${bam}"}" TMP/tmp.collate |\\
	samtools fastq \\
		-n \\
		--threads 1 \\
		-1 TMP/${prefix}.paired.R1.fq.gz \\
		-2 TMP/${prefix}.paired.R2.fq.gz \\
		-s TMP/${prefix}.singleton.fq.gz \\
		-0 TMP/${prefix}.other.fq.gz



mv -v TMP/${prefix}*.fq.gz ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
stub:
"""
touch ${meta.id}.paired.R1.fq.gz
touch ${meta.id}.paired.R2.fq.gz
touch ${meta.id}.singleton.fq.gz
touch ${meta.id}.other.fq.gz
touch versions.yml
"""
}

