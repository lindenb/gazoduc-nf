include {moduleLoad} from '../../modules/utils/functions.nf'

def TAG="VISTA_ENHANCER"

workflow ANNOTATE_VISTA {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		source_ch =  DOWNLOAD(genomeId)
		annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)
	emit:
		output = annotate_ch.output
		count = annotate_ch.count
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
	def url = genome.vista_enhancers_url
        def reference = genome.fasta
        def whatis = "Vista Enhancers from ${url}"
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit bedtools")}

set -o pipefail
wget -O - "${url}" |\
        gunzip -c |\
        cut -f2-5 |\
        java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\
        bgzip > TMP/${TAG}.bed.gz

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > ${TAG}.header
"""
}

process ANNOTATE {
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	path(tabix)
	path(tbi)
	path(header)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}" --merge-logic '${TAG}:unique' -O b -o TMP/${TAG}.bcf '${vcf}'
bcftools index TMP/${TAG}.bcf

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
