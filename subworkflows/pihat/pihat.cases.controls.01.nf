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
nextflow.enable.dsl=2


include {moduleLoad} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {PIHAT01} from './pihat.01.nf'
include	{VCF_INTER_CASES_CONTROLS_01} from '../bcftools/vcf.inter.cases.controls.01.nf'

workflow PIHAT_CASES_CONTROLS_01 {
	take:
		genomeId
		vcf
		cases
		controls
	main:
		version_ch = Channel.empty();
		
		inter_ch = VCF_INTER_CASES_CONTROLS_01([:], vcf, cases, controls)
		version_ch = version_ch.mix(inter_ch.version)

		pihat_ch = PIHAT01(genomeId,vcf,inter_ch.all_samples)
		version_ch = version_ch.mix(pihat_ch.version)


		remove_ch= REMOVE_SAMPLES(pihat_ch.pihat_removed_samples ,cases,controls)
		version_ch = version_ch.mix(remove_ch.version)

		version_ch = MERGE_VERSION("pihat",version_ch.collect())

	emit:
		version = version_ch
		cases = remove_ch.cases
		controls = remove_ch.controls
		pihat_png = pihat_ch.pihat_png
		removed_samples = pihat_ch.pihat_removed_samples
		plink_genome = pihat_ch.plink_genome
		pihat_sample2avg_png  = pihat_ch.pihat_sample2avg_png
	}

process REMOVE_SAMPLES {
executor "local"
input:
	val(removed_samples)
	val(cases)
	val(controls)
output:
	path("filtered.cases.txt"),emit:cases
	path("filtered.controls.txt"),emit:controls
	path("version.xml"),emit:version
script:
"""
hostname 1>&2

cut -f 1 "${removed_samples}" | sort -T . > a
sort -T . "${cases}" > b
sort -T . "${controls}" > c


comm -13 a b > filtered.cases.txt
comm -13 a c > filtered.controls.txt

test -s filtered.cases.txt
test -s filtered.controls.txt


#############################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">remove sample with high pihat</entry>
	<entry key="removed samples">\$(cat a | paste -s -d ' ')</entry>
	<entry key="count(cases)">\$(wc -l < filtered.cases.txt)</entry>
	<entry key="count(controls)">\$(wc -l < filtered.controls.txt)</entry>
</properties>
EOF


rm a b c
sleep 10
"""
}
