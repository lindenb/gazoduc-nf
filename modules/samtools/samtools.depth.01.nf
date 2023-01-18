/*

Copyright (c) 2023 Pierre Lindenbaum

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

include {assertNotEmpty;getKeyValue;moduleLoad;isBlank} from '../../modules/utils/functions.nf'

process SAMTOOLS_DEPTH_01 {
tag "${row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.sample}.depth.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:""
	def bed = row.bed==null || row.bed.name.equals("NO_FILE")?"":"-M -L \"${row.bed}\" "
	def ref = row.reference?"--reference \"${row.reference}\" ":""
	def has_bed = !isBlank(bed)
	def mapq = row.mapq?:1

"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

${has_bed?"samtools view -F 3844 --min-MQ \"${mapq}\" --uncompressed -O BAM ${bed} ${ref} \"${row.bam}\" |\\":""}
samtools depth -a -q "${mapq}" ${has_bed?" -b \"${row.bed}\" ":""} ${isBlank(row.interval)?"-":" -r \"${row.interval}\" \"${row.bam}\""}  |\
awk -F '\t' 'BEGIN{N=0;T=0.0;} {T+=int(\$3);N++;} END {printf("#sample\tbam\treference\tinterval\tbed\tcoverage\\n");printf("${row.sample}\t${row.bam}\t${row.reference?:""}\t${row.interval?:"."}\t${row.bed?:"."}\t%f\\n",(N==0?0.0:T/N));}' > "${row.sample}.depth.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools depth</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="interval">${row.interval?:"."}</entry>
	<entry key="bed">${row.bed?:"."}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}
