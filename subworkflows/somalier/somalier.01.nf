include {isHg19;isHg38;hasFeature;getModules} from './functions.nf'
include {vcfSamples} from './vcfsamples.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow SOMALIER_VCF_01 {
	take:
		meta
		reference
		vcf
		pedigree
	input:
		version_ch = Channel.empty()

		exe_ch = DOWNLOAD_SOMALIER(meta)
		version_ch = version_ch.mix(somalier_ch.version)

		sites_ch = DOWNLOAD_SITES(meta,reference)
		version_ch = version_ch.mix(sites_ch.version)
	
		vcf_samples_ch = vcfSamples(meta,vcf,pedigree)	
		version_ch = version_ch.mix(vcf_samples_ch.version)
	
		somalier_ch = APPLY_SOMALIER(meta, exe_ch.executable ,sites_ch.vcf , vcf_samples_ch.samples,vcf_samples_ch.pedigree, vcf, reference)
		version_ch = version_ch.mix(somalier_ch.version)

		version_ch = MERGE_VERSION(meta, "somalier", "somalier", version_ch.collect())
	emit:
		version = version_ch
	}


process APPLY_SOMALIER {
afterScript "rm -rf extracted TMP"
input:
	val(meta)
	val(somalier)
	val(sites)
        path(samples)
	path(pedigree)
	val(vcf)
	val(reference)
output:
	path("*vcf.somalier.archive.zip"),emit:zip
script:
def prefix = meta.prefix?:"SOMALIER."

if(hasFeature("somalier"))
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
set -x

mkdir TMP

bcftools query -f '%CHROM\t%POS0\t%END\\n' "${sites}" > TMP/jeter.bed

if [ ! -z "${vcf.endsWith(".list")?"Y":""}" ] ; then

	bcftools concat --file-list "${vcf}" --regions-file TMP/jeter.bed -O u  --allow-overlaps --remove-duplicates |\
		bcftools view --samples-file "${samples}" -O b -o TMP/jeter.bcf

else

	bcftools view --regions-file TMP/jeter.bed --samples-file "${samples}" -O b -o TMP/jeter.bcf

fi

bcftools index TMP/jeter.bcf

# remove header
tail -n +2 "${pedigree}" > TMP/jeter.ped
test -s TMP/jeter.ped

${somalier} extract -d extracted \
	--sites "${sites}" \
	-f "${reference}" TMP/jeter.bcf

find \${PWD}/extracted -type f -name "*.somalier" | xargs ${somalier} relate --output-prefix=${prefix}vcf -p TMP/jeter.ped

# may not exist
touch "${prefix}vcf.groups.tsv"
zip -9 "${prefix}vcf.somalier.archive.zip" ${prefix}vcf.groups.tsv ${prefix}vcf.html ${prefix}vcf.pairs.tsv ${prefix}vcf.samples.tsv

rm -rf extracted
"""
}
