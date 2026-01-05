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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {UTR_ANNOTATOR_DOWNLOAD_01} from '../../modules/vep/utrannotator.download.01.nf'
def TAG="VEP"

workflow ANNOTATE_VEP {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

	
		if(params.genomes[genomeId].containsKey("vep_module") && params.genomes[genomeId].containsKey("ucsc_name") && (params.genomes[genomeId].ucsc_name.equals("hg38") || params.genomes[genomeId].ucsc_name.equals("hg19"))) {
			source_ch = UTR_ANNOTATOR_DOWNLOAD_01([:],genomeId)
			annotate_ch = ANNOTATE(genomeId, source_ch.output,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			}
		else
			{
			out1 = vcfs
			out2 = Channel.empty()
			}
	emit:
		output = out1
		count = out2
	}



process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(vep_utr)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	def genome = params.genomes[genomeId]
	def vep = vep_utr.toRealPath()
"""
hostname 1>&2
${moduleLoad("bcftools")}
module load ${genome.vep_module}
mkdir -p TMP

# vep

set +u
export PERL5LIB=\${PERL5LIB}:`dirname ${vep}`

bcftools view "${vcf}" |\
	${genome.vep_invocation.replace("--cache", "--plugin UTRannotator,${vep} --cache") } |\
	bcftools sort -T TMP/x -O b -o TMP/${TAG}.bcf

bcftools index --force TMP/${TAG}.bcf

###
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
