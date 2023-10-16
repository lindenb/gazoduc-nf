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

include {moduleLoad;getKeyValue;hasFeature} from '../../modules/utils/functions.nf'
include {VALIDATE_CASE_CONTROL_PED_01} from '../../modules/pedigree/validate.case.ctrl.pedigree.01.nf'
include {VCF_INTER_CASES_CONTROLS_01} from '../bcftools/vcf.inter.cases.controls.01.nf'
include {PEDIGREE_FOR_RVTESTS} from '../../modules/rvtests/rvtests.cases.controls.ped.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {PIHAT_CASES_CONTROLS_01} from '../pihat/pihat.cases.controls.01.nf'

workflow BURDEN_SAMPLES_PART_01 {
	take:
		meta
		reference
		vcf
		pedigree
	main:
		to_zip = Channel.empty()
		version_ch = Channel.empty()
		
		ped_ch = VALIDATE_CASE_CONTROL_PED_01(meta,pedigree)
		version_ch = version_ch.mix(ped_ch.version)

		vcf_inter_ch = VCF_INTER_CASES_CONTROLS_01(meta.plus(["with_tabix":true]),vcf,
					ped_ch.cases_list, ped_ch.controls_list )
		version_ch = version_ch.mix(vcf_inter_ch.version)

		if(hasFeature(meta,"pihat")) {
			pihat = PIHAT_CASES_CONTROLS_01(meta,reference,file(vcf),ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(pihat.version)
			to_zip = to_zip.mix(pihat.pihat_png)
			to_zip = to_zip.mix(pihat.pihat_sample2avg_png)
			to_zip = to_zip.mix(pihat.removed_samples)
			to_zip = to_zip.mix(pihat.plink_genome)

			rebuild_ch = REBUILD_PEDIGREE(meta,pedigree,pihat.cases,pihat.controls)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}
		else
			{
			rebuild_ch = REBUILD_PEDIGREE(meta,pedigree,ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}

                rvtests_ped_ch = PEDIGREE_FOR_RVTESTS(meta,new_ped_ch)
                version_ch = version_ch.mix(rvtests_ped_ch.version)

		version_ch = MERGE_VERSION(meta, "PedigreeBurden", "Pedigree for burden", version_ch.collect())

	emit:
		version = version_ch
		zip = to_zip
		pedigree = new_ped_ch
		rvtest_pedigree = rvtests_ped_ch.pedigree
	}


process REBUILD_PEDIGREE {
executor "local"
input:
	val(meta)
	path(pedigree)
	path(cases)
	path(controls)
output:
	path("case_control.ped"),emit:pedigree
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail

sort -T . -t '\t' -k2,2 "${pedigree}" > tmp1.ped
cat  "${cases}" "${controls}" | sort | uniq > tmp2.ped

join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6' tmp1.ped  tmp2.ped |\
	sort -T . -t '\t' -k1,1 -k2,2  > case_control.ped
test -s case_control.ped
rm tmp1.ped tmp2.ped



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">make pedigree after samples have been removed.</entry>
	<entry key="pedigree">${pedigree.toRealPath()}</entry>
	<entry key="cases">${cases.toRealPath()}</entry>
	<entry key="controls">${controls.toRealPath()}</entry>
	<entry key="samples.removed">\$(comm -23 <(cut -f 2 "${pedigree}" | sort)   <(cut -f 2 "case_control.ped" | sort)  | paste -s -d ' ' )</entry>
</properties>
EOF
"""
}
