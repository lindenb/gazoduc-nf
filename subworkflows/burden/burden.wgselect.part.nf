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

include {moduleLoad;getKeyValue;hasFeature} from '../../modules/utils/functions.nf'
include {BED_CLUSTER_01} from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {WGSELECT_01} from '../wgselect/wgselect.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow BURDEN_WGSELECT_PART_01 {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:
		to_zip = Channel.empty()
		version_ch = Channel.empty()	

                cluster_ch = BED_CLUSTER_01(meta.plus(
                        ["bed_cluster_method":getKeyValue(meta,"bed_cluster_method","--size 5mb")]),
                        reference,
                        bed
                        )
		version_ch = version_ch.mix(cluster_ch.version)

		wgselect_ch = WGSELECT_01(meta, reference, vcf, pedigree, cluster_ch.output.splitText().map{it.trim()}.map{T->file(T)})
		version_ch = version_ch.mix(wgselect_ch.version.first())
		to_zip = to_zip.mix(wgselect_ch.variants_list)

		
		version_ch = MERGE_VERSION(meta, "burdenWgSelect", "Burden WgSelectPart", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

	emit:
		version = version_ch
		zip = to_zip
                variants_list = wgselect_ch.variants_list /** file containing count of variants at each step */
                contig_vcfs = wgselect_ch.contig_vcfs /** path to all vcf concatenated per contigs */
                vcfs =  wgselect_ch.vcfs /** path to all chunks of vcf */

	}
