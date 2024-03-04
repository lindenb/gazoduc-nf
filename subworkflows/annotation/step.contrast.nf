include {moduleLoad;isBlank} from '../../modules/utils/functions.nf'
def TAG="CONTRAST"

workflow ANNOTATE_ELSWHERE {
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
	def args=""
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP

# sample in VCF
bcftools query -l "${vcf}" | sort -T TMP | uniq > TMP/samples.txt

# samples in ped with phenotype
awk '{if(length(\$6) >0) print \$2}' '${pedigree}' | sort | uniq > TMP/jeter2.txt

# common samples
comm -12 TMP/samples.txt TMP/jeter2.txt > TMP/common.txt

# new pedigree
tr -s " " "\t" < '${pedigree}' | tr -s "\t" |\\
	sort -T TMP -t '\t' -k2,2 |\\
	join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6' - TMP/common.txt > TMP/jeter.ped


awk -F '\t' '{if(\$6=="case" || \$6=="2") print \$2;}' TMP/jeter.ped > TMP/cases.txt
awk -F '\t' '{if(\$6=="control" || \$6=="1") print \$2;}' TMP/jeter.ped > TMP/controls.txt

if test -s TMP/cases.txt && test -s TMP/controls.txt
then
	bcftools +contrast -0 TMP/controls.txt -1 TMP/cases.txt ${args} -O b -o  TMP/${TAG}.bcf 
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
