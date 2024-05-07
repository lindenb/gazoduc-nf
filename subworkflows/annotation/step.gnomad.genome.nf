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
include {hasFeature;isSoftFilter;isBlank;backDelete} from './annot.functions.nf'
def TAG="GNOMAD_GENOME"

workflow ANNOTATE_GNOMAD_GENOME {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:

	
		if(hasFeature("gnomad_genome") && !isBlank(params.genomes[genomeId],"gnomad_genome") ) {
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
   	def url = genome.clinvar_vcf_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Annotation with gnomad genome <code>${genome.gnomad_genome}</code>. max AF :${params.annotations.gnomad_max_AF}, Population: <b>${params.annotations.gnomad_population}</b></dd>
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
	def max_AF = params.annotations.gnomad_max_AF
	def pop =  params.annotations.gnomad_population

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP OUTPUT
set -o pipefail

bcftools view "${row.vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \
			--bufferSize 10000 \
			--max-af ${max_AF} \
			--gnomad "${genome.gnomad_genome}" --fields "${pop}" |\
	bcftools view -O b -o TMP/${TAG}.bcf


 	if  ${!isSoftFilter("gnomad")} ; then
        bcftools view  -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_RF"' -O b -o TMP/jeter.bcf TMP/${TAG}.bcf
        mv TMP/jeter.bcf TMP/${TAG}.bcf
    fi


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
