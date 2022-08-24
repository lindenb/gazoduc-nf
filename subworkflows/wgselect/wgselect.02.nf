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
include {getBoolean;getKeyValue;getModules;getGnomadExomePath;getGnomadGenomePath;isHg19;isHg38;hasFeature;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_PER_CONTIG_01} from '../bcftools/bcftools.concat.contigs.01.nf'
include {VCF_TO_BED_01} from '../../modules/jvarkit/jvarkit.vcf2bed.01.nf'
include {TO_LIST_01} from  '../../modules/utils/to_list.nf'
include {WGSELECT_01} from './wgselect.01.nf'

workflow WGSELECT_02 {
	take:
		meta
		reference
		vcf
		pedigree
	main:
		version_ch = Channel.empty()
		list_ch = TO_LIST_01(meta,vcf)
		version_ch = version_ch.mix(list_ch.version)
		
		vcfs_ch = list_ch.output.splitText().map{it.trim()}
		
		tobed_ch = VCF_TO_BED_01(meta,vcfs_ch)
		version_ch = version_ch.mix(tobed_ch.version)

		each_bed = tobed_ch.output.splitCsv(header:true,sep:'\t').map{T->T.bed}
	
		wch1_ch = WGSELECT_01(meta, reference, vcf, pedigree, each_bed)
		version_ch = version_ch.mix(wch1_ch.version)
		
		version_ch = MERGE_VERSION(meta, "WGSelect", "WG Select", version_ch.collect())

	emit:
		version = version_ch /** version */
                variants_list = wch1_ch.variants_list /** file containing count of variants at each step */
               	contig_vcfs = wch1_ch.contig_vcfs /** path to all vcf concatenated per contigs */
                vcfs = wch1_ch.vcfs /** path to all chunks of vcf */


	}

