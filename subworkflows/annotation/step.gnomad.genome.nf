include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
def TAG="GNOMAD_GENOME"

workflow ANNOTATE_GNOMAD_GENOME {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

	
		if(params.genomes[genomeId].containsKey("gnomad_genome") ) {
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
	def max_AF = 0.01
	def pop =  "AF_NFE"

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP OUTPUT


bcftools view "${row.vcf}" |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad \
			--bufferSize 10000 \
			--max-af ${max_AF} \
			--gnomad "${genome.gnomad_genome}" --fields "${pop}" |\
	bcftools view -O b -o TMP/${TAG}.bcf

bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "\${PWD}/OUTPUT/${TAG}.bed"
}
EOF

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv TMP/${TAG}.* OUTPUT/
rm -fv "${row.vcf}" "${row.index}"
"""
}
