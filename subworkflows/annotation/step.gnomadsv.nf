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
include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../gnomad/download_gnomad_sv.01.nf'


def TAG="GNOMADSV"

workflow ANNOTATE_GNOMADSV {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:
		if(hasFeature("gnomadsv")) {
			source_ch =  DOWNLOAD_GNOMAD_SV_01([:], genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.index, vcfs)

			out1 = annotate_ch.output
			out2 = annotate_ch.count
			out3 = MAKE_DOC(genomeId).output
			}
		else
			{
			out1 = vcfs
			out2 = Channel.empty()
			out3 = Channel.empty()
			}
	emit:
		output = out1
		count = out2
		doc = out3
}


process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def genome = params.genomes[genomeId]
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Frequent SV in Gnomad. min-AF:${params.annotations.gnomadsv_min_AF} Field ${params.annotations.gnomadsv_population} </dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(tabix)
	path(tbi)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def AF=params.annotations.gnomadsv_min_AF?:0.1
	def pop = params.annotations.gnomadsv_population?:"POPMAX_AF"
	def whatis = "GNOMAD SV with ${pop} > ${AF}"
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools htslib")}
mkdir -p TMP OUTPUT

# check header exists
gunzip -c "${tabix}" | head -n 1 | tr "\t" "\\n" | grep -F '${pop}'

set -o pipefail

tabix --print-header --region "${row.bed}" "${tabix}" |\
	awk '(NR==1) {C=-1;for(i=1;i<=NF;i++) if(\$i=="${pop}") C=i;next;} {if(\$C!="NA" && \$C*1.0 > ${AF} ) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4);}'  |\
        LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
        bgzip > TMP/gnomad.sv.bed.gz

tabix  --force -p bed -f TMP/gnomad.sv.bed.gz

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > TMP/header.txt


bcftools annotate -a "${tabix}" -h TMP/header.txt -c "CHROM,FROM,TO,${TAG}"  --merge-logic '${TAG}:unique'  -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT
${backDelete(row)}
"""
}
