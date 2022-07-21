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

include {isBlank;moduleLoad} from '../../modules/utils/functions.nf'

process GRAPHTYPER_GENOTYPE_01 {
tag "${row.interval}"
cpus 1
afterScript "rm -rf results TMP"
input:
	val(meta)
	val(graphtyper)
	val(row)
output:
	tuple val(row),path("genotyped.bcf"),emit:output
	path("version.xml"),emit:version
script:
	if(!row.reference) {exit 1,"reference missing"}
	if(!row.bams) {exit 1,"bams missing"}
	if(!row.interval) {exit 1,"interval missing"}

	def avg_cov_by_readlen= row.avg_cov_by_readlen?:""
	def arg2 = isBlank(avg_cov_by_readlen)?"":"--avg_cov_by_readlen=${avg_cov_by_readlen}"
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir TMP

export TMPDIR=\${PWD}/TMP

${graphtyper} genotype \
	"${row.reference}" \
	--force_no_copy_reference \
	--force_use_input_ref_for_cram_reading \
	--sams=${row.bams} \
	--region=${row.interval} \
	--threads=${task.cpus} \
	${arg2}

find \${PWD}/results/ -type f -name "*.vcf.gz" | grep -v '/input_sites/' > TMP/vcf.list

bcftools concat --file-list TMP/vcf.list \
	--allow-overlaps --remove-duplicates \
	--threads ${task.cpus} -O b -o "genotyped.bcf"

bcftools index --threads ${task.cpus} "genotyped.bcf"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run graphtyper</entry>
</properties>
EOF
"""
}
