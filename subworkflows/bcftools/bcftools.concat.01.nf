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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {BCFTOOL_CONCAT_FILE_LIST_01} from '../../modules/bcftools/bcftools.concat.file.list.01.nf'
include {BCFTOOL_CONCAT_COLLECT_01} from '../../modules/bcftools/bcftools.concat.collect.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow BCFTOOLS_CONCAT_01 {
take:
	meta
	vcfs
main:
	version_ch = Channel.empty()
	d1_ch = SQRT_FILE(meta, vcfs)

	d2_ch = d1_ch.clusters.splitText().map{it.trim()}

	d3_ch = BCFTOOL_CONCAT_FILE_LIST_01(meta, d2_ch)
	version_ch = version_ch.mix(d3_ch.version)

	d4_ch = BCFTOOL_CONCAT_COLLECT_01(meta, d3_ch.vcf.collect() )
	version_ch = version_ch.mix(d4_ch.version)

	 version_ch = MERGE_VERSION(meta, "concat", "concat vcfs", version_ch.collect())
emit:
	vcf = d4_ch.vcf
	version = version_ch
}


