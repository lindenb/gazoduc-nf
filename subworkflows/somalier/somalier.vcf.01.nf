include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
//include {vcfSamples} from './vcfsamples.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../modules/somalier/somalier.download.sites.nf'
include {VCF_INTER_PED_01} from '../../modules/bcftools/vcf.inter.pedigree.01.nf'



workflow SOMALIER_VCF_01 {
	take:
		meta
		reference
		vcf
		pedigree
	main:
		assertFileExists(reference,"reference must be defined")
		assertFileExists(vcf,"vcf must be defined")
		version_ch = Channel.empty()

		exe_ch = DOWNLOAD_SOMALIER(meta)
		version_ch = version_ch.mix(exe_ch.version)

		sites_ch = SOMALIER_DOWNLOAD_SITES(meta,reference)
		version_ch = version_ch.mix(sites_ch.version)
	
		if(!isBlank(pedigree)) {
			vcf_samples_ch = VCF_INTER_PED_01(meta.plus(["pedigree_type","other"]),vcf,pedigree)	
			version_ch = version_ch.mix(vcf_samples_ch.version)
			ped_ch = vcf_samples_ch.pedigree
			}
		else
			{
			ped_ch = ""
			}
	
		somalier_ch = APPLY_SOMALIER(meta, reference, exe_ch.executable ,sites_ch.vcf ,vcf, ped_ch)
		version_ch = version_ch.mix(somalier_ch.version)

		version_ch = MERGE_VERSION(meta, "somalier", "somalier on vcf", version_ch.collect())
	emit:
		version = version_ch
		zip = somalier_ch.zip
	}


process APPLY_SOMALIER {
afterScript "rm -rf extracted TMP"
input:
	val(meta)
	val(reference)
	val(somalier)
	val(sites)
	val(vcf)
	val(pedigree)
output:
	path("${prefix}somalier.vcf.zip"),emit:zip
	path("version.xml"),emit:version
script:
	prefix = meta.prefix?:""

"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
set -x

mkdir TMP

bcftools query -f '%CHROM\t%POS0\t%END\\n' "${sites}" > TMP/jeter.bed

## extract samples
if [ ! -z "${pedigree}" ] ; then
	# remove header
	grep -v "^fid\tiid"  "${pedigree}" > TMP/jeter.ped
	test -s TMP/jeter.ped

	cut -f 2 TMP/jeter.ped | sort | uniq > TMP/jeter.samples.txt
	test -s TMP/jeter.samples.txt
fi


if [ ! -z "${vcf.endsWith(".list")?"Y":""}" ] && [ -f TMP/jeter.samples.txt ] ; then

	bcftools concat --file-list "${vcf}" --regions-file TMP/jeter.bed -O u  --allow-overlaps --remove-duplicates |\
		bcftools view --samples-file TMP/jeter.samples.txt -O b -o TMP/jeter.bcf

else if [ ! -z "${vcf.endsWith(".list")?"Y":""}" ] ; then

	bcftools concat --file-list "${vcf}" --regions-file TMP/jeter.bed -O b  --allow-overlaps --remove-duplicates -o TMP/jeter.bcf

else if [ -f TMP/jeter.samples.txt ] ; then

	bcftools view --regions-file TMP/jeter.bed --samples-file "${samples}" -O b -o TMP/jeter.bcf

else

	bcftools view --regions-file TMP/jeter.bed -O b -o TMP/jeter.bcf

fi

bcftools index TMP/jeter.bcf


${somalier} extract -d extracted \
	--sites "${sites}" \
	-f "${reference}" TMP/jeter.bcf

mkdir "${prefix}somalier.vcf"

find \${PWD}/extracted -type f -name "*.somalier" | xargs ${somalier} relate --output-prefix=${prefix}somalier.vcf/${prefix}vcf ${pedigree.isEmpty()?"":"-p TMP/jeter.ped"}

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run somalier on vcf</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="pedigree">${pedigree}</entry>
</properties>
EOF



# may not exist
touch "${prefix}somalier.vcf/${prefix}vcf.groups.tsv"
zip -9 -r "${prefix}somalier.vcf.zip" "${prefix}somalier.vcf"

"""
}
