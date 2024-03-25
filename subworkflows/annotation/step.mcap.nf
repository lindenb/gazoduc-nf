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
include {hasFeature;isBlank;backDelete;isHg19;isHg18} from './annot.functions.nf'

def TAG="MCAP"
def WHATIZ="M-CAP is the first pathogenicity classifier for rare missense variants in the human genome that is tuned to the high sensitivity required in the clinic. By combining previous pathogenicity scores (including SIFT, Polyphen-2 and CADD) with novel features and a powerful model, we attain the best classifier at all thresholds, reducing a typical exome/genome rare (<1%) missense variant (VUS) list from 300 to 120, while never mistaking 95% of known pathogenic variants as benign."

String getUrl(genomeId) {
	if(isHg19(genomeId)) return "http://bejerano.stanford.edu/mcap/downloads/dat/mcap_v1_4.txt.gz"
	return "";
	}

workflow ANNOTATE_MCAP {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

             if(hasFeature("mcap") && !getUrl(genomeId).isEmpty()) {
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
memory "2g"
input:
	val(genomeId)
output:
	path("${TAG}.bed.gz"),emit:bed
	path("${TAG}.bed.gz.tbi"),emit:tbi
	path("${TAG}.header"),emit:header
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = getUrl(genomeId)
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit")}

wget -O TMP/jeter.tsv.gz "${url}"

gunzip -c TMP/jeter.tsv.gz |\
	tail -n +2 |\
	cut -f 1-5 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
		bgzip > TMP/${TAG}.bed.gz && \
	tabix -p bed -f TMP/${TAG}.bed.gz


mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=1,Type=Float,Description="${WHATIZ} ${url}">' > ${TAG}.header
"""
}



process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def url = getUrl(genomeId)
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>${WHATIZ} <a href="${url}">${url}</a></dd>
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

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,POS,REF,ALT,${TAG}" -O b --merge-logic '${TAG}:max' -o TMP/${TAG}.bcf '${row.vcf}'
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
