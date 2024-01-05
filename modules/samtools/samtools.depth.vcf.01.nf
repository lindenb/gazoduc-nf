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

include {assertNotEmpty;getKeyValue;moduleLoad;isBlank} from '../../modules/utils/functions.nf'

process SAMTOOLS_DEPTH_VCF_01 {
tag "${row.bam}"
memory "3g"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("coverage.bcf"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:""
	def bed = row.bed==null || row.bed.name.equals("NO_FILE")?"":"-M -L \"${row.bed}\" "
	def ref = row.reference?"--reference \"${row.reference}\" ":""
	def has_bed = !isBlank(bed)
	def mapq = row.mapq?:1
	def bam = row.bam
"""
hostname 1>&2
${moduleLoad("samtools bcftools")}
set -o pipefail


	mkdir TMP
	samtools samples "${bam}" | awk -F '\t' '{printf("#CHROM\tPOS\t%s\\n",\$1);}' > TMP/jeter.tsv

	${has_bed?"samtools view -F 3844 --min-MQ \"${mapq}\" --uncompressed -O BAM ${bed} ${ref} \"${bam}\" |\\":""}
	samtools depth -a -q "${mapq}" ${has_bed?" -b \"${row.bed}\" ":""} ${isBlank(row.interval)?"-":" -r \"${row.interval}\" \"${bam}\""}  >> TMP/jeter.tsv

	# convert to vcf
	awk -F '\t' 'BEGIN {printf("##fileformat=VCFv4.2\\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\\"depth\\">\\n##INFO=<ID=MINDP,Number=1,Type=Integer,Description=\\"min depth\\">\\n##INFO=<ID=MAXDP,Number=1,Type=Integer,Description=\\"max depth\\">\\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\\"depth\\">\\n");} /^#CHROM/{printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\\n",\$3);next;} {printf("%s\t%s\t.\tN\t.\t.\t.\tDP=%d;MINDP=%d;MAXDP=%d\tDP\t%s\\n",\$1,\$2,\$3,\$3,\$3,\$3);}' TMP/jeter.tsv > TMP/jeter.vcf
	bcftools view ${task.cpus?"--threads ${task.cpus}":""} -o TMP/jeter.vcf.gz -O z TMP/jeter.vcf

	bcftools reheader --threads ${task.cpus} -T TMP --fai "${row.reference}.fai" -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz

	bcftools sort --max-mem "${task.memory.giga}G" -T TMP -O b -o TMP/jeter2.bcf TMP/jeter2.vcf.gz
	bcftools index --threads ${task.cpus} TMP/jeter2.bcf

	bcftools view --threads ${task.cpus} ${has_bed?"--regions-file \"${row.bed}\"":""}  ${isBlank(row.interval)?"":"--regions \"${row.interval}\""} -O b -o coverage.bcf TMP/jeter2.bcf
	bcftools index --threads ${task.cpus} coverage.bcf


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools depth en convert to VCF</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="interval">${row.interval?:"."}</entry>
	<entry key="bed">${row.bed?:"."}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}
