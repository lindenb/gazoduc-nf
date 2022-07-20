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

include {assertNotEmpty;getKeyValue;moduleLoad} from '../../modules/utils/functions.nf'

process SAMTOOLS_FASTQ_01 {
tag "${file(row.bam).name} ${row.contig?:""} ${row.interval?:""} ${row.bed?:""}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	val(row)
output:
	path("output.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def bam = assertNotEmpty(row.bam,"bam must be defined")
	def reference = assertNotEmpty(row.reference,"reference must be defined")
	def sample = assertNotEmpty(row.sample,"sample must be defined")
	def bed = row.bed?"--region-file \"${row.bed}\"":""
	def contig = row.contig?"\"${row.contig}\"":""
	def interval = row.interval?"\"${row.interval}\"":""
	def has_region = !bed.isEmpty() || !contig.isEmpty()  ||  !interval.isEmpty()
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail
mkdir TMP

# get unmapped reads from input, convert to fastq
${has_region?"samtools view --reference \"${reference}\" -M --uncompressed ${bed} \"${bam}\" ${contig} ${interval} |\\":""}
	samtools collate -f --threads 4 -O -u --no-PG --reference "${reference}" ${has_region?"-":"\"${bam}\""} TMP/tmp.collate |\
	samtools fastq -N --threads 1 -1 TMP/jeter.R1.fq.gz -2 TMP/jeter.R2.fq.gz -s TMP/jeter.Rx.fq.gz -0 TMP/jeter.R0.fq.gz -n

gunzip -c  TMP/jeter.Rx.fq.gz | paste - - - - | awk -F '\t' '(\$1 ~ /\\/1\$/) {OFS="\t";gsub(/\\/1\$/,"",\$1);print}' | gzip --best > TMP/unpaired.R1.tsv.gz
gunzip -c  TMP/jeter.Rx.fq.gz | paste - - - - | awk -F '\t' '(\$1 ~ /\\/2\$/) {OFS="\t";gsub(/\\/2\$/,"",\$1);print}' | gzip --best > TMP/unpaired.R2.tsv.gz



mv TMP/jeter.R1.fq.gz "${sample}.R1.fq.gz"
mv TMP/jeter.R2.fq.gz "${sample}.R2.fq.gz"
mv TMP/jeter.R0.fq.gz "${sample}.R0.fq.gz"
mv TMP/unpaired.R1.tsv.gz "${sample}.unpaired.R1.tsv.gz"
mv TMP/unpaired.R2.tsv.gz "${sample}.unpaired.R2.tsv.gz"


echo "sample\tbam\treference\tR1\tR2\tother\tunpairedR1\tunpairedR2" > output.tsv
echo "${sample}\t${bam}\t${reference}\t\${PWD}/${sample}.R1.fq.gz\t\${PWD}/${sample}.R2.fq.gz\t\${PWD}/${sample}.R0.fq.gz\t\${PWD}/${sample}.unpaired.R1.tsv.gz\t\${PWD}/${sample}.unpaired.R2.tsv.gz" >> output.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${bam}</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
	<entry key="bed">${bed}</entry>
	<entry key="contig">${contig}</entry>
	<entry key="interval">${interval}</entry>
</properties>
EOF
"""
stub:
"""
touch sample2bam.tsv
echo "<properties/>" > version.xml
"""
}

