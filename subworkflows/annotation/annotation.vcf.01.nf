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
nextflow.enable.dsl=2


include {moduleLoad;getKeyValue} from '../../modules/utils/functions.nf'

workflow ANNOT_VCF_01 {
	take:
		meta
		reference
		vcfs
	main:
		annot_ch = ANNOT(meta,reference,vcfs)
	emit:
		vcf = annot_ch.vcf

}

process ANNOT {
tag "${file(vcf).name}"
afterScript "rm -f jeter.bcf"
memory "5g"
input:
	val(meta)
	val(reference)
	val(vcf)
output:
	path("vep.bcf"),emit:vcf
	path("vep.bcf.csi"),emit:csi
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}

module load ensembl-vep/104.3

bcftools view "${vcf}" |\
vep --cache  --format vcf --force_overwrite --output_file STDOUT --no_stats --offline \
	--dir_cache /LAB-DATA/BiRD/resources/apps/vep  --species homo_sapiens --cache_version 104 \
	--assembly GRCh37 --merged --fasta "${reference}" --use_given_ref --vcf |\
bcftools view -O b -o vep.bcf
bcftools index vep.bcf


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">annot with VEP</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}
