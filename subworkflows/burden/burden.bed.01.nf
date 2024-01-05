/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {BED_CLUSTER_01} from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from '../wgselect/wgselect.01.nf'
include {RVTESTS01_VCF_01} from '../../modules/rvtests/rvtests.vcf.01.nf'
include {RVTESTS_POST_PROCESS} from '../rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {PIHAT_CASES_CONTROLS_01} from '../pihat/pihat.cases.controls.01.nf'
include {BURDEN_SAMPLES_PART_01} from './burden.samples.part.nf'

workflow BURDEN_BED_01 {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:
		to_zip = Channel.empty()
		version_ch = Channel.empty()
	
		samples_ch = BURDEN_SAMPLES_PART_01(meta,reference,vcf,pedigree)
		version_ch = version_ch.mix(samples_ch.version)
		to_zip  = to_zip.mix(samples_ch.zip)


                cluster_ch = BED_CLUSTER_01(meta.plus(
                        ["bed_cluster_method":getKeyValue(meta,"bed_cluster_method","--size 5mb")]),
                        reference,
                        bed
                        )
		version_ch = version_ch.mix(cluster_ch.version)

		wgselect_ch = WGSELECT_01(meta, reference, vcf, samples_ch.pedigree, cluster_ch.output.splitText().map{it.trim()}.map{T->file(T)})
		version_ch = version_ch.mix(wgselect_ch.version.first())
		to_zip = to_zip.mix(wgselect_ch.variants_list)


		assoc_ch = RVTESTS01_VCF_01(meta, reference,wgselect_ch.vcfs.splitText().map{it.trim()}, samples_ch.rvtest_pedigree )
		version_ch = version_ch.mix(assoc_ch.version.first())

		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference,file(vcf).name,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden in BED", "Burden in BED.", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

	emit:
		version = version_ch
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
