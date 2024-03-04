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
def TAG="ALPHAMISSENSE"

workflow ANNOTATE_ALPHAMISSENSE {
	take:
		genomeId
		vcfs /** json: vcf,vcf_index */
	main:

	
		if(hasFeature("alphamissense") && !isBlank(params.genomes[genomeId],"alphamissense_url")) {
			source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
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


process DOWNLOAD {
afterScript "rm -rf TMP"
memory "3g"
input:
        val(genomeId)
output:
	path("${TAG}.tsv.gz"),emit:bed
	path("${TAG}.tsv.gz.tbi"),emit:tbi
	path("${TAG}.header"),emit:header
script:
	def genome = params.genomes[genomeId]
   	def url = genome.alphamissense_url
    def whatis="Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492"

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit htslib")}
set -o pipefail
mkdir -p TMP

set -o pipefail
mkdir -p TMP
set -x

wget -O TMP/jeter.tsv.gz "${url}"

gunzip -c TMP/jeter.tsv.gz |\
	grep -v '^#' |\
	cut -f 1 | uniq | sort -T TMP | uniq |\
	awk '{printf("%s\t%s\\n",\$1,\$1);}' |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${genome.reference}" --column 2 --convert SKIP |\
	awk -F '\t' '{printf("s|^%s\t|%s\t|\\n",\$1,\$2);}' > TMP/jeter.sed


gunzip -c TMP/jeter.tsv.gz |\
	grep -v '^#'  |\
	cut -f 1-4,9,10 |\
	sed -f TMP/jeter.sed |\
	LC_ALL=C sort --buffer-size=1000M -T TMP -t '\t' -k1,1 -k2,2n |\
	uniq |\
	bgzip >  TMP/${TAG}.tsv.gz && \

tabix -s 1 -b 2 -e 2  TMP/${TAG}.tsv.gz

mv TMP/${TAG}.tsv.gz ./
mv TMP/${TAG}.tsv.gz.tbi ./

echo '##INFO=<ID=${TAG}_PATHOGENOCITY,Number=1,Type=Float,Description="${whatis}. ${url}.">' >  ${TAG}.header
echo '##INFO=<ID=${TAG}_CLASS,Number=.,Type=String,Description="${whatis}. ${url}.">' >> ${TAG}.header

"""
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
<dd>Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492 </dd>
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
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,POS,REF,ALT,${TAG}_PATHOGENOCITY,${TAG}_CLASS" --merge-logic "${TAG}_PATHOGENOCITY:max,{TAG}_CLASS:unique" -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
