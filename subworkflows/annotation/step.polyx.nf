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

def TAG="POLYX"

workflow ANNOTATE_POLYX {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

		 if(hasFeature("polyx") && (params.annotations.polyx as int)>0) {
			    annotate_ch = ANNOTATE(genomeId,vcfs)
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
	def url = genome.rmsk_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Poly-X with jvarkit </dd>
</dl>
__EOF__
"""
}


process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(genomeId)
	//tuple path(vcf),path(vcf_idx),path(bed)
	path(json)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def genome = params.genomes[genomeId]
	def row = slurpJsonFile(json)

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP OUTPUT
set -o pipefail

bcftools view "${row.vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfpolyx \
			-n '${params.annotations.polyx}' \
			--reference '${genome.fasta}' \
			--tag POLYX   |\
	bcftools view -O b -o TMP/${TAG}.bcf

bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
