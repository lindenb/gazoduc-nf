/*

Copyright (c) 2022 Pierre Lindenbaum

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


include {moduleLoad;getKeyValue;getGnomadExomePath;getGnomadGenomePath} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {PIHAT01} from './pihat.01.nf'
include	{VCF_INTER_CASES_CONTROLS_01} from '../bcftools/vcf.inter.cases.controls.01.nf'

workflow PIHAT_CASES_CONTROLS_01 {
	take:
		meta
		reference
		vcf
		cases
		controls
	main:
		version_ch = Channel.empty();
		
		inter_ch = VCF_INTER_CASES_CONTROLS_01(meta, vcf, cases, controls)
		version_ch = version_ch.mix(inter_ch.version)

		pihat_ch = PIHAT01(meta,reference,vcf,inter_ch.all_samples)
		version_ch = version_ch.mix(pihat_ch.version)

		version_ch = MERGE_VERSION(meta, "pihat", "pihat on ${vcf} with cases and controls.", version_ch.collect())

	emit:
		version = version_ch
		/* pihat_png = pihat_ch.pihat_png
		removed_samples = pihat_ch.pihat_removed_samples
		plink_genome = pihat_ch.plink_genome
		pihat_pdf  = pihat_ch.pihat_pdf */
	}

/*
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">prepare pihat data per sexual contig</entry>
	<entry key="vcf">${vcf}</entry>
	<entry key="contig">${contig}</entry>
</properties>
EOF
*/