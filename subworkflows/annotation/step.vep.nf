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
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
