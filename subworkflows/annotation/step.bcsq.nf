include {moduleLoad} from '../../modules/utils/functions.nf'

def TAG="BCSQ"

workflow ANNOTATE_BCSQ {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		if(params.genomes[genomeId].containsKey("gff3")) {
			annotate_ch = ANNOTATE(genomeId,vcfs)
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
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	def genome= params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

                
bcftools csq -O b --force --local-csq --ncsq 10000 --fasta-ref "${reference}" --gff-annot "${genome.gff3}" -o TMP/${TAG}.bcf '${vcf}'
bcftools index TMP/${TAG}.bcf

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
