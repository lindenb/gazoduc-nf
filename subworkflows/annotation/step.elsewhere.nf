include {moduleLoad;isBlank} from '../../modules/utils/functions.nf'
def TAG="ELSWHERE"

workflow ANNOTATE_ELSWHERE {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

		if(!isBlank(params.other_vcfs)) {
			annotate_ch = ANNOTATE(genomeId, file(params.other_vcfs),vcfs)
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
tag "${vcf.name} ${bed.name}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(others)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP


echo -n '##INFO=<ID=${TAG},Number=.,Type=String,Description="Samples found in other files: '>  TMP/${TAG}.header
sort "${others}" | paste -sd ' ' | tr -d '\\n' >>   TMP/${TAG}.header
echo '">' >> TMP/${TAG}.header


bcftools query -l "${vcf}" | sort -T TMP | uniq > TMP/samples.1.txt


cat "${others}" | while read F
do
	# samples in other file
	bcftools query -l "\${F}" |  sort -T TMP | uniq > TMP/samples.2.txt
	
	# remove common samples
	comm -13 TMP/samples.1.txt TMP/samples.2.txt > TMP/samples.3.txt
	if test -s TMP/samples.3.txt
	then
		bcftools view --regions-file "${vcf}" --samples-file TMP/samples.3.txt -O u "\${F}" |\\
		bcftools view -c 1 -O u |\\
		bcftools norm -f '${reference}' --multiallelics -any -O u |\\
		bcftools query -i 'ALT!="*"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\\n]' |\
		awk -F '\t' '{G=\$NF;if(G ~ /^[0\\.][|/][0\\.]\$/ || G=="0" || G==".") next; print;}' |\
		cut -f 1-5 |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4 > TMP/jeter2.bed

		# merge with previous bed if any
		if test -s TMP/jeter.bed
		then
			LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4  --merge TMP/jeter.bed TMP/jeter2.bed | uniq > TMP/jeter3.bed
			mv TMP/jeter3.bed TMP/jeter.bed
			rm TMP/jeter2.bed
		else
			mv TMP/jeter2.bed TMP/jeter.bed
		fi	
	done
done

if test -s TMP/jeter.bed
then
	bgzip TMP/jeter.bed
	tabix --force -s 1 -b 2 -e 2 TMP/jeter.bed.gz
	bcftools annotate -a TMP/jeter.bed.gz --columns 'CHROM,POS,REF,ALT,${TAG}' --header-lines TMP/${TAG}.header --merge-logic '${TAG}:unique' -O b -o  TMP/${TAG}.bcf 
else

	bcftools view -O b -o TMP/${TAG}.bcf "${vcf}"

fi

bcftools index TMP/${TAG}.bcf

mv  TMP/${TAG}.bcf ./
mv  TMP/${TAG}.bcf.csi ./

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/count.tsv
mv TMP OUTPUT
"""
}
