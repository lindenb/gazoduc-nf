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

include {getVersionCmd;jvarkit;isHg19;isHg38;isBlank;hasFeature;moduleLoad;getGnomadGenomePath;getGnomadExomePath} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {ANNOTATE_VCF_01} from './annotation.vcf.01.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../subworkflows/jvarkit/jvarkit.vcf2intervals.nf'

workflow ANNOTATE_VCF_02 {
	take:
		meta
		reference
		vcf
		bed /* limit to that BED or NO_FILE */
	main:
		version_ch = Channel.empty()

		rch = JVARKIT_VCF_TO_INTERVALS_01(meta, reference, vcf, bed )
		version_ch = version_ch.mix(rch.version)

		interval_ch = rch.bed.splitCsv(header:false,sep:'\t').
			map{T->[interval:T[0]+":"+((T[1] as int)+1)+"-"+T[2],vcf:T[3]]}
	
		annotate_ch = ANNOTATE_VCF_01(meta , reference, interval_ch)
		version_ch = version_ch.mix(annotate_ch.version)

	
		file_list_ch = COLLECT_TO_FILE_01([:],annotate_ch.bedvcf.
				splitCsv(header:false,sep:'\t').
				map{T->T[1]}.collect())
	
		concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output)
		version_ch = version_ch.mix(concat_ch.version)

	
		version_ch = MERGE_VERSION(meta, "annot vcf", "annotation of VCF", version_ch.collect())


	emit:
		version = version_ch
	}

