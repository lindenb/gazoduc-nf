include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'

/* ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA).  */

def TAG="REMAP"

workflow ANNOTATE_REMAP {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		if(params.genomes[genomeId].containsKey("remap_url")) {
			source_ch =  DOWNLOAD(genomeId)
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
		output = out1
		count = out2
}

process DOWNLOAD{
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
	def url = genome.remap_url
        def reference = genome.fasta
        def whatis = "ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA). ${url}"
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
        gunzip -c |\
        cut -f1,2,3 |\
        java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\
	bedtools merge |\
	sed 's/\$/\t1/' |\
        bgzip > TMP/${TAG}.bed.gz

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
	//tuple path(vcf),path(vcf_idx),path(bed)
	path(json)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	def row = slurpJsonFile(json)	
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}"  -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT

rm -fv "${row.vcf}" "${row.index}"
"""
}
