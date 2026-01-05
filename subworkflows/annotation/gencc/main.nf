/*

Copyright (c) 2026 Pierre Lindenbaum

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


 
workflow GENCC {
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

		ANNOTATE(DOWNLOAD.out.output, vcfs)
		versions = versions.mix(ANNOTATE.out.versions)
		
	emit:
		vcf = ANNOTATE.out.vcf
		versions
}

process DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:header
	path("versions.yml"),emit:versions
script:
	def TAG="GENCC"
	def URL = "https://search.thegencc.org/download/action/submissions-export-tsv"
	def WHATIZ = "gencc The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly https://thegencc.org/about.html"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP



cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz
1:${k1.hg19}\t${base}/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv
join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
    sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed --columns "gtf.feature,gene_name" -R "${fasta}" "${gtf}" |\\
	awk -F '\t' '(\$4=="gene")' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T TMP -t '\t' -k4,4 |\\
	uniq > TMP/genes.bed

wget --no-check-certificate -q  -O - "${URL}" |\\
	awk -F '\t' '(NR==1) {C=-1;D=-1;for(i=1;i<=NF;i++) {if(\$i=="\\"gene_symbol\\"") C=i;if(\$i=="\\"disease_title\\"") D=i;}next;} {if(C<1 ||  D<1) next; G=\$C;H=\$D; gsub(/[^A-Za-z0-9\\.\\-]+/,"_",G);gsub(/[^A-Za-z0-9\\.\\-]+/,"_",H);  printf("%s\t%s\\n",G,H);}'  |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1  | uniq > TMP/org.tsv


LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' TMP/genes.bed TMP/org.tsv |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > TMP/${TAG}.bed.gz


tabix --force -p bed TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ} ${URL}.">' > ${TAG}.header

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${URL}"
END_VERSIONS
"""
}


process ANNOTATE {
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(tabix),path(tbi),path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
	tuple path("*.bcf"),path("*.bcf.csi"),emit:bed
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:vcf.baseName+".gencc"
	def TAG="GENCC"
"""
mkdir -p TMP OUTPUT

bcftools annotate \\
	--threads ${task.cpus} \\
	-a "${tabix}" \\
	-h "${header}" \\
	-c "CHROM,FROM,TO,${TAG}" \\
	--merge-logic '${TAG}:unique' \\
	-O b \\
	-o TMP/${prefix}.bcf '${vcf}'

mv TMP/${prefix}.bcf ./
mv TMP/${prefix}.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
