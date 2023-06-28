/*

Copyright (c) 2023 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include {VCF_INTER_SAMPLES_01 as INTER_CASES}   from '../../modules/bcftools/vcf.inter.samples.01.nf' addParams(prefix : "cases")
include {VCF_INTER_SAMPLES_01 as INTER_CONTROLS} from '../../modules/bcftools/vcf.inter.samples.01.nf' addParams(prefix : "controls")
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'

/**

find common samples in a case|control file and a pdf.
Check there is no overlap
check no file is empty 

*/
workflow VCF_INTER_CASES_CONTROLS_01 {
	take:
		vcf
		cases_file
		controls_file
	main:
		version_ch = Channel.empty()
		cases_ch = INTER_CASES(vcf,cases_file)
		version_ch = version_ch.mix(cases_ch.version)

		ctrls_ch = INTER_CONTROLS(vcf,controls_file)
		version_ch = version_ch.mix(ctrls_ch.version)
		

		
		
		uniq_ch = UNIQ(
				cases_ch.common,
				cases_ch.vcf_only,
				ctrls_ch.common,
				ctrls_ch.vcf_only
				)
		version_ch = version_ch.mix(uniq_ch.version)

		version_ch = MERGE_VERSION("intersection samples in VCF cases and controls", version_ch.collect())
	emit:
		version = version_ch
		cases = cases_ch.common
		controls = ctrls_ch.common
		vcf_only = uniq_ch.vcf_only
		all_samples = uniq_ch.all_samples
	}


process UNIQ {
tag "${common_cases} ${common_controls}"
executor "local"
afterScript "rm -f a.txt b.txt c.txt"
input:
        val(common_cases)
       	val(vcfonly_cases)
        val(common_controls)
	val(vcfonly_controls)
output:
	path("all_samples.txt"),emit:all_samples
	path("vcf_only.txt"),emit:vcf_only
	path("version.xml"),emit:version
script:
"""
hostname 1>12
set -o pipefail

# all samples
cat "${common_cases}" "${common_controls}" | sort -T . | uniq > all_samples.txt

# check no common samples between cases and controls
sort "${common_cases}" > a.txt
sort "${common_controls}" > b.txt
comm -12 a.txt b.txt > c.txt

echo -n "Common samples cases/controls (should be empty):" 1>&2
paste -s -d ' ' < c.txt 1>&2
test ! -s c.txt
rm c.txt 

# in vcf, absent from both files case/controls
cat "${vcfonly_cases}" "${vcfonly_controls}" | sort -T . | uniq -d > vcf_only.txt

################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="id">${task.process}</entry>
	<entry key="description">merge case/controls</entry>
</properties>
EOF
"""
}	

