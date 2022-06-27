include {VCF_INTER_SAMPLES_01 as INTER_CASES; VCF_INTER_SAMPLES_01 as INTER_CONTROLS} from '../../modules/bcftools/vcf.inter.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow VCF_INTER_CASES_CONTROLS_01 {
	take:
		meta
		vcf
		cases_file
		controls_file
	main:
		version_ch = Channel.empty()
		cases_ch = INTER_CASES(meta,vcf,cases_file)
		version_ch = version_ch.mix(cases_ch.version)

		version_ch = Channel.empty()
		ctrls_ch = INTER_CONTROLS(meta,vcf,controls_file)
		version_ch = version_ch.mix(ctrls_ch.version)
		
		uniq_ch = UNIQ(meta,
				cases_ch.common,
				cases_ch.vcf_only,
				ctrls_ch.common,
				ctrls_ch.vcf_only
				)
		version_ch = version_ch.mix(uniq_ch.version)

		version_ch = MERGE_VERSION(meta, "vcf inter cases/ctrls", "intersection samples in VCF cases and controls", version_ch.collect())
	emit:
		version = version_ch.version
		cases = uniq_ch.cases
		controls = uniq_ch.controls
		vcf_only = uniq_ch.vcf_only
		all_samples = uniq_ch.all_samples
	}


process UNIQ {
executor "local"
input:
	val(meta)
	path(common_cases)
	path(vcfonly_cases)
	path(common_controls)
	path(vcfonly_controls)
output:
	path(common_cases),emit:cases
	path(common_controls),emit:controls
	path("vcf_only.txt"),emit:vcf_only
	path("all_samples.txt"),emit:all_samples
	path("version.xml"),emit:version
script:
"""

# all samples
cat "${common_cases}" "${common_controls}" | sort -T . | uniq > all_samples

# check no common samples between cases and controls
comm -12 <(sort "${common_cases}") <( sort "${common_controls}") > c.txt

echo -n "Common samples cases/controls (should be empty):" 1>&2
cat c.txt 1>&2
test ! -s c.txt
rm c.txt

# in vcf, absent from both files case/controls
cat "${vcfonly_cases}" "${vcfonly_controls}" | sort -T . | uniq -d > vcf_only.txt

cat << EOF
<properties id="${task.process}">
	<entry key="id">${task.process}</entry>
	<entry key="description">merge case/controls</entry>
</properties>
EOF

"""
}
