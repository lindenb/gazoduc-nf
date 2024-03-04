include {moduleLoad;isBlank} from '../../modules/utils/functions.nf'
def TAG="DNM2"

workflow ANNOTATE_DNM2 {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

		if(!isBlank(params.pedigree)) {
			annotate_ch = ANNOTATE(genomeId, file(params.pedigree),vcfs)
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
	path(pedigree)
	tuple path(vcf),path(vcf_idx),path(bed)
output:
	tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/count.tsv"),emit:count
script:
	def args="--use-NAIVE "
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

# sample in VCF
# new pedigree
tr -s " " "\t" < '${pedigree}' | tr -s "\t" |\\
	cut -f1-4 |\
	sort -T TMP -t '\t' -k2,2 > TMP/tmp.ped

# get samples with parent
awk -F '\t' '(\$3!="0" && \$4!="0")' TMP/tmp.ped | cut -f2,3,4 | tr "\t" "\n" | sort -T TMP | uniq > TMP/jeter.txt

# common samples
join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4' TMP/tmp.ped TMP/jeter.txt  > TMP/tmp2.ped
mv TMP/tmp2.ped TMP/tmp.ped


if test -s TMP/tmp.ped
then
	bcftools +trio-dnm2 -P TMP/tmp.ped ${args} -O b -o  TMP/${TAG}.bcf 
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
