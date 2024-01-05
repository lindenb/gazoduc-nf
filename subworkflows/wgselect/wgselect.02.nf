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

include {moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from './wgselect.01.nf'
include {LINUX_SPLIT} from '../../modules/utils/split.nf'
include {JVARKIT_GATK_HARD_FILTERING_01} from '../jvarkit/jvarkit.gatk_hard_filtering.01.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../jvarkit/jvarkit.vcf2intervals.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'

workflow WGSELECT_02 {
	take:
		genomeId
		vcf
		pedigree
		bed /* limit to that bed */
	main:
		version_ch = Channel.empty()
		
		if( (params.wgselect.gatk_hardfiltering_percentile as double) > 0 ) {
			vcf2bed_ch = VCF_TO_BED([with_header: false],vcf)
			version_ch = version_ch.mix(vcf2bed_ch.version)
		
			in_ch = vcf2bed_ch.bed.splitCsv(header:false,sep:'\t').map{T->[
					interval: T[0]+":"+((T[1] as int)+1)+"-"+T[2],
					vcf: T[3]
					]}
		
			hard_filters_ch = JVARKIT_GATK_HARD_FILTERING_01([percentile: params.wgselect.gatk_hardfiltering_percentile], in_ch)
			
			version_ch = version_ch.mix(hard_filters_ch.version)
			
			hard_filters = hard_filters_ch.output
			}
		else
			{
			hard_filters = file("NO_FILE")
			}
		
		tobed_ch = JVARKIT_VCF_TO_INTERVALS_01([with_header:false,distance: params.wgselect.distance ,min_distance: params.wgselect.min_distance], vcf, bed)
		version_ch = version_ch.mix(tobed_ch.version)

		wch1_ch = WGSELECT_01(
			[:],
			genomeId,
			tobed_ch.bed.splitCsv(header:false,sep:'\t').
				map{T->[
					interval: T[0]+":"+((T[1] as int)+1)+":"+T[2],
					vcf: T[3],
					pedigree: pedigree,
					]}.
				combine(hard_filters).
				map{T->T[0].plus(hard_filters:T[1])}
				
			)

		version_ch = version_ch.mix(wch1_ch.version)		
		version_ch = MERGE_VERSION("WGSelect", version_ch.collect())

	emit:
		version = version_ch /** version */
                variants_list = wch1_ch.variants_list /** file containing count of variants at each step */
               	contig_vcfs = wch1_ch.contig_vcfs /** path to all vcf concatenated per contigs */
                vcfs = wch1_ch.vcfs /** path to all chunks of vcf */


	}


