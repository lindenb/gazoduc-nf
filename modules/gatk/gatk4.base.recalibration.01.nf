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
include {moduleLoad;isBlank} from './../utils/functions.nf'

process GATK4_BASE_RECALIBRATOR_01 {
cache 'lenient'
memory '5g'
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	val(row)
output:
        tuple val(row),path("${row.sample}.recal.table"),emit:output
        path("version.xml"),emit:version
script:
	def ks0 = (row.recalibration_vcfs?:"").split("[, ]+").
			findAll{T->!T.isEmpty()}.collect{V->" --known-sites "+V}.
			join(" ")
	def ks = isBlank(ks0) && !isBlank(meta.dbsnp) ? "--known-sites \"${meta.dbsnp}\"":ks0
	def bed = (row.bed?"-L \"${row.bed}\"":"")
"""
hostname 1>&2
${moduleLoad("gatk4 bcftools")}

mkdir TMP

# when testing, provide an empty VCF
if [ ! -z "${isBlank(ks)?"Y":""}" ] ; then

cat << EOF | bcftools view -O z -o TMP/jeter.vcf.gz
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF

bcftools index --tbi TMP/jeter.vcf.gz
fi

# allele specific annotation are not supported in non-gvcf mode
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" BaseRecalibrator \
	-I "${row.bam}" \
	 ${bed} \
	 ${isBlank(ks)?"--known-sites TMP/jeter.vcf.gz":ks} \
	-O "${row.sample}.recal.table" \
	-R "${row.reference}"

###########################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">gatk BaseRecalibrator</entry>
	<entry key="bams">${row.bam}</entry>
	<entry key="bed">${bed}</entry>
	<entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF
"""

stub:
"""
touch "${row.sample}.recal.table"
echo "<properties/>" > version.xml
"""
}
