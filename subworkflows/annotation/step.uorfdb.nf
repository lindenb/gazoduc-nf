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
include {hasFeature;isBlank;backDelete;isHg38} from './annot.functions.nf'
def TAG="UORFDB"
def WHATIZ="upstream open reading frame database (uORF)"
workflow ANNOTATE_UORFDB {
	take:
		genomeId
		bed
		vcfs /** json : tuple vcf,vcf_index */
	main:

	
		if(hasFeature("uorfdb") && isHg38(genomeId) ) {
			source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			out3 = MAKE_DOC().output
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


process DOWNLOAD {
afterScript "rm -rf TMP"
memory "3g"
input:
        val(genomeId)
output:
	path("${TAG}.bed.gz"),emit:bed
	path("${TAG}.bed.gz.tbi"),emit:tbi
	path("${TAG}.header"),emit:header
script:
	def hg="hg38"
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
       	def url = "https://www.bioinformatics.uni-muenster.de/tools/uorfdb/download/uORF_dump_uORFdb.tsv"
"""
hostname 1>&2
${moduleLoad("htslib jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\\
	awk -F '\t' '(\$2=="${hg}" && \$22!="" && \$23!="" && \$30!="")' |\\
	cut -f 3,22,23,30 |\\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\\
	LC_ALL=C sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip 	> TMP/${TAG}.bed.gz
	

tabix -p bed -f TMP/${TAG}.bed.gz

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="Kozak_strength uORF ${WHATIZ} ${url}">' > ${TAG}.header

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./
"""
}


process MAKE_DOC {
executor "local"
output:
	path("${TAG}.html"),emit:output
script:
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>${WHATIZ}</dd>
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
	path(header)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}"  --merge-logic '${TAG}:unique' -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index --force TMP/${TAG}.bcf

bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

mv -v TMP/${TAG}.* OUTPUT
${backDelete(row)}
"""
}
