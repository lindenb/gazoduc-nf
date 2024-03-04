include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {DOWNLOAD_GNOMAD_SV_01} from '../gnomad/download_gnomad_sv.01.nf'


def TAG="GNOMADSV"

workflow ANNOTATE_GNOMADSV {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		if(params.genomes[genomeId].containsKey("remap_url")) {
			source_ch =  DOWNLOAD_GNOMAD_SV_01([:], genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.index, vcfs)

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
input:
	path(tabix)
	path(tbi)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def AF=0.1
	def pop = "POPMAX_AF"
	def whatis = "GNOMAD SV with ${pop} > ${AF}"
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools htslib")}
mkdir -p TMP OUTPUT

tabix --print-header --region "${row.bed}" "${tabix}" |\
	awk '(NR==1) {C=-1;for(i=1;i<=NF;i++) if(\$i=="${pop}") C=i;next;} {if(\$C!="NA" && \$C*1.0 > ${AF} ) printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4);}'  |\
        LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
        bgzip > TMP/gnomad.sv.bed.gz

tabix  --force -p bed -f TMP/gnomad.sv.bed.gz

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > TMP/header.txt


bcftools annotate -a "${tabix}" -h TMP/header.txt -c "CHROM,FROM,TO,${TAG}"  --merge-logic '${TAG}:unique'  -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "\${PWD}/OUTPUT/${TAG}.bed"
}
EOF

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT
rm -fv "${row.vcf}" "${row.index}"
"""
}
