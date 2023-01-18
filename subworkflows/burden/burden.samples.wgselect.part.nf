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
include {BED_CLUSTER_01} from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from '../wgselect/wgselect.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {PIHAT_CASES_CONTROLS_01} from '../pihat/pihat.cases.controls.01.nf'
include {BURDEN_SAMPLES_PART_01} from './burden.samples.part.nf'
include {BURDEN_WGSELECT_PART_01} from './burden.wgselect.part.nf'

/**
 *
 * select samples AND apply wgselect for the new pedigree
 *
 */
workflow BURDEN_SAMPLE_WGSELECT_PART_01 {
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

		wgsel_ch = BURDEN_WGSELECT_PART_01(meta, reference, vcf, samples_ch.pedigree, bed)
		to_zip  = to_zip.mix(wgsel_ch.zip)
		version_ch = version_ch.mix(wgsel_ch.version)
		
		version_ch = MERGE_VERSION(meta, "burdenSamplesWgSelect", "Burden Samples and WgSelect Part", version_ch.collect())
	emit:
		version = version_ch
		zip = to_zip
                pedigree = samples_ch.pedigree
                rvtest_pedigree = samples_ch.rvtest_pedigree 

                variants_list = wgsel_ch.variants_list /** file containing count of variants at each step */
                contig_vcfs = wgsel_ch.contig_vcfs /** path to all vcf concatenated per contigs */
                vcfs =  wgsel_ch.vcfs /** path to all chunks of vcf */


	}

