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

include {parseBoolean;isBlank;moduleLoad} from '../../modules/utils/functions.nf'

process BWA_MEM_01 {
tag "${row.sample} ${row.R1} ${row.R2?:""}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.sample}.sorted.bam"), path("${row.sample}.sorted.bam.bai"),emit:output
	path("version.xml"),emit:version
script:
	if(!row.containsKey("R1")) throw new IllegalArgumentException("R1 missing");
	if(!row.containsKey("sample")) throw new IllegalArgumentException("sample missing");
	if(!row.containsKey("genomeId")) throw new IllegalArgumentException("genomeId missing");
	def sample = row.sample
	def genome = params.genomes[row.genomeId]
	def reference = genome.fasta
	def R1 = row.R1
	def R2 = row.R2?:""
	def is_interleaved = parseBoolean(row.interleaved?:"false")
	def do_index = parseBoolean(row.with_index?:"true")
	def ID = row.ID?:sample
	def CN = row.CN?:"Nantes"
	def PL = row.PL?:"ILLUMINA"
	def LB = row.LB?:sample
	def cpus1 = ((task.cpus?:2) as int)
	def cpus2 = cpus1 -1
"""
hostname 1>&2
${moduleLoad("samtools bwa")}
set -o pipefail
mkdir -p TMP
set -x

bwa mem ${is_interleaved?"-p":""} -t ${cpus2} -R '@RG\\tID:${ID}\\tSM:${sample}\\tLB:${LB}.R0\\tCN:${CN}\\tPL:${PL}' "${reference}" "${R1}" ${isBlank(R2)?"":"\"${R2}\""} |\
	samtools view -O BAM -o TMP/jeter.bam

if ${!isBlank(R2)} ; then

	# collate
	samtools collate -l 5 --threads ${cpus1}  --output-fmt BAM -u --no-PG --reference "${reference}" -T TMP/tmp.collate -o TMP/jeter2.bam TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam

	# fixmate
	samtools fixmate --threads  ${cpus1} -mc --output-fmt BAM TMP/jeter.bam TMP/jeter2.bam
	mv TMP/jeter2.bam TMP/jeter.bam

fi


samtools sort --threads ${cpus1} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
mv TMP/jeter2.bam TMP/jeter.bam

samtools index  -@ ${cpus1} TMP/jeter.bam

mv "TMP/jeter.bam" "${sample}.sorted.bam"
mv "TMP/jeter.bam.bai" "${sample}.sorted.bam.bai"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">map fastq to reference with bwa</entry>
	<entry key="sample">${sample}</entry>
	<entry key="reference">${reference}</entry>
	<entry key="interleaved">${is_interleaved}</entry>
	<entry key="R1">${R1}</entry>
	<entry key="R2">${R2}</entry>
	<entry key="ID">${ID}</entry>
	<entry key="CN">${CN}</entry>
	<entry key="PL">${PL}</entry>
	<entry key="bwa.version">\$(bwa 2>&1 | grep Version)</entry>
</properties>
EOF
"""
}

