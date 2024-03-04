include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'

def TAG="RMSK"

workflow ANNOTATE_RMSK {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

             	if(params.genomes[genomeId].containsKey("rmsk_url")) {
                        source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
                        out1 = annotate_ch.output
                        out2 = annotate_ch.count
                        }
                else
                    	{
                        out1 = vcfs
                        out2 = Channel.empty()
                        }
	emit:
		output = annotate_ch.output
		count = annotate_ch.count
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
	def url = genome.rmsk_url
	def whatis="Repeat Masker from ${url}"
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
	gunzip -c |\
	cut -f6-8 |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
		sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		sed 's/\$/\t1/' |\
		bgzip > TMP/${TAG}.bed.gz && \
	tabix -p bed -f TMP/${TAG}.bed.gz


mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="${whatis}">' > ${TAG}.header
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
	path("OUTPUT/count.tsv"),emit:count
script:
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}" -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf

rm -fv "${row.vcf}" "${row.index}"


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT

"""
}
