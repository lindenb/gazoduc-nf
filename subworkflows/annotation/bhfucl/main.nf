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



workflow ANNOTATE_BHFUCL {
	take:
		meta
		fasta
		fai
		dict
		gtf
		vcfs /* meta, vcf,vcf_index */
	main:
		versions = Channel.empty()

		DOWNLOAD(fasta,fai,dict,gtf)
		versions = versions.mix(DOWNLOAD.out.versions)

		ANNOTATE(
			DOWNLOAD.out.bed_in,
			DOWNLOAD.out.bed_out,
			vcfs
			)
		versions = versions.mix(ANNOTATE.out.versions)
	emit:
		vcf = ANNOTATE.out.vcf
		versions
	}
		
process DOWNLOAD {
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf),path(gtf_tbi)
	path("versions.yml"),emit:versions
output:
	tuple val(meta1),path("*IN.bed.gz"), path("*IN.bed.gz.tbi"), path("*IN.header"), emit:bed_in
	tuple val(meta1),path("*OUT.bed.gz"),path("*OUT.bed.gz.tbi"),path("*OUT.header"),emit:bed_out
script:
    def TAG = task.ext.tag?:"BHFUCL"
    def URL = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
    def WHATIZ = "Cardiovascular Gene Ontology Annotation Initiative ${URL}"
	def extend = task.ext.extend?:1000
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

wget -O - "${URL}" |\\
	gunzip -c |\\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\\
	uniq | LC_ALL=C sort -T TMP | uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}_IN.bed.gz

tabix --force -p bed TMP/${TAG}_IN.bed.gz

gunzip -c TMP/${TAG}_IN.bed.gz |\\
	bedtools slop -b ${extend} -g ${fai} |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	uniq |\\
	bgzip  > TMP/${TAG}_OUT.bed.gz

tabix --force -p bed TMP/${TAG}_OUT.bed.gz


mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ} ${url}">' > ${TAG}_IN.header
echo '##INFO=<ID=${TAG}_NEAR,Number=.,Type=String,Description="Near gene distance=${extend}. ${WHATIZ} ${url}">' > ${TAG}_OUT.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${URL}"
END_VERSIONS
"""
}



process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(bed_in),path(tabix_in),path(header_in)
	tuple val(meta2),path(bed_out),path(tabix_out),path(header_out)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:output
script:
    def TAG = task.ext.tag?:"BHFUCL"
	def prefix=task.ext.prefix?:vcf.baseName+".bhfucl"
    def distance = task.ext.distance?:1000;
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP OUTPUT



bcftools annotate \\
	--threads ${task.cpus} \\
	-a "${bed_in}" \\
	-h "${header_in}" \\
	--write-index \\
	-c "CHROM,POS,END,${TAG}" \\
	-O b \\
	--merge-logic '${TAG}:unique' \\
	-o TMP/jeter.bcf \\
	'${vcf}'

bcftools annotate \\
	--threads ${task.cpus} \\
	-a "${tabix_out}" \\
	-h "${header_out}" \\
	--write-index \\
	--keep-sites -e 'INFO/${TAG} != ""' \\
	-c "CHROM,POS,END,${TAG}_NEAR" \\
	-O b \\
	--merge-logic '${TAG}_NEAR:unique' \\
	-o TMP/${prefix}.bcf \\
	TMP/jeter.bcf



mv TMP/*${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
