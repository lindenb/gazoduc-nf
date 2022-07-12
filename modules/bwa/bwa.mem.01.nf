/*

Copyright (c) 2022 Pierre Lindenbaum

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

include {parseBoolean;assertNotEmpty;getKeyValue;moduleLoad} from '../../modules/utils/functions.nf'

process BWA_MEM_01 {
tag "${row.sample}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	val(row)
output:
	tuple val("${row.sample}"),path("${row.sample}.sorted.bam"),emit:bam
	path("${row.sample}.sorted.bam.bai"),emit:bai
	path("version.xml"),emit:version
script:
	def sample = row.sample
	def reference = row.reference
	def R1 = row.R1
	def R2 = row.R2?:""
	def is_interleaved = parseBoolean(row.interleaved?:"false")
	def do_index = parseBoolean(row.with_index?:"true")
	def ID = row.ID?:sample
	def CN = row.CN?:"Nantes"
	def PL = row.PL?:"ILLUMINA"
	def LB = row.LB?:sample
"""
hostname 2>&1
${moduleLoad("samtools bwa")}
set -o pipefail
mkdir TMP

bwa mem ${is_interleaved?"-p":""} -t ${task.cpus} -R '@RG\\tID:${ID}\\tSM:${sample}\\tLB:${LB}.R0\\tCN:${CN}\\tPL:${PL}' "${reference}" "${R1}" ${R2.isEmpty()?"\"${R2}\"":""} |\
	samtools view -O BAM -o TMP/jeter.bam


samtools sort --threads ${task.cpus} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
mv TMP/jeter2.bam TMP/jeter.bam

samtools index  -@ ${task.cpus} TMP/jeter.bam

mv "TMP/jeter.bam" "${sample}.sorted.bam"
mv "TMP/jeter.bam.bai" "${sample}.sorted.bam.bai"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="R1">${R1}</entry>
	<entry key="R2">${R2}</entry>
</properties>
EOF
"""
stub:
"""
touch sample2bam.tsv
echo "<properties/>" > version.xml
"""
}

