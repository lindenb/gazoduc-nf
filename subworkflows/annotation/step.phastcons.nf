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

def TAG="PHASTCONS"
def WHATIZ="PhastCons is a hidden Markov model-based method that estimates the probability that each nucleotide belongs to a conserved element, based on the multiple alignment. It considers not just each individual alignment column, but also its flanking columns."

String getUrl(genomeId) {
	if(isHg19(genomeId)) return "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons46way/vertebrate.phastCons46way.bw"
	if(isHg38(genomeId)) return "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.phastCons470way.bw"
	return "";
	}
workflow ANNOTATE_HMC {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

             if(hasFeature("hmc") && (isHg19(genomeId) || isHg38(genomeId))) {
                        source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.output,vcfs)
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
	path("${TAG}.bw"),emit:output
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url = getUrl(genomeId)
"""
hostname 1>&2
mkdir -p TMP

wget -O TMP/${TAG}.bw "${url}"
mv -v TMP/${TAG}.bw ./
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
memory "3g"
input:
	path(bw)
	path(json)
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP OUTPUT


bcftools view '${row.vcf}' |\\
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfbigwig -B '${bw}' --tag ${TAG} |\\
bcftools view -O b -o TMP/${TAG}.bcf

bcftools index TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
