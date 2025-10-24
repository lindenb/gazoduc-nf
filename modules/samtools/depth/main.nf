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

include {assertNotEmpty;getKeyValue;isBlank} from '../../modules/utils/functions.nf'

process SAMTOOLS_DEPTH{
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
cpus 1
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(optional_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.depth.tsv"),emit:depth
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:bam.baseName
	def bed = row.bed==null || row.bed.name.equals("NO_FILE")?"":"-M -L \"${row.bed}\" "
	def ref = fasta?"--reference \"${fasta}\" ":""
	def has_bed = optional_bed?true:false
	def mapq = task.ext.mapq?:1

"""
hostname 1>&2
set -o pipefail

${has_bed?"samtools view -F 3844 --min-MQ \"${mapq}\" --uncompressed -O BAM ${optional_bed} ${ref} \"${row.bam}\" |\\":""}
samtools depth -a -q "${mapq}" ${has_bed?" -b \"${row.bed}\" ":""} ${isBlank(row.interval)?"-":" -r \"${row.interval}\" \"${row.bam}\""}  |\
awk -F '\t' 'BEGIN{N=0;T=0.0;} {T+=int(\$3);N++;} END {printf("#sample\tbam\treference\tinterval\tbed\tcoverage\\n");printf("${row.sample}\t${row.bam}\t${row.reference?:""}\t${row.interval?:"."}\t${row.bed?:"."}\t%f\\n",(N==0?0.0:T/N));}' > "${row.sample}.depth.tsv"

cat << END_VERSIONS > versions.yml
"${task.process}":
	samtools: "\$(samtools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
