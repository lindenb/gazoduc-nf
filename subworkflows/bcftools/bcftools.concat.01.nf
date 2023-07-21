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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {BCFTOOL_CONCAT_FILE_LIST_01 as CONCAT1} from '../../modules/bcftools/bcftools.concat.file.list.01.nf' addParams(prefix:"tmp.")
include {BCFTOOL_CONCAT_FILE_LIST_01 as CONCAT2} from '../../modules/bcftools/bcftools.concat.file.list.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


workflow BCFTOOLS_CONCAT_01 {
take:
	meta
	vcfs /* path containing the path to a set of indexed vcf files */
	bed /* or NO_FILE */
main:
	version_ch = Channel.empty()
	d1_ch = SQRT_FILE(min_file_split:meta.min_file_split?:20, suffix:".list", vcfs)

	d2_ch = d1_ch.output.splitText().map{file(it.trim())}

	d3_ch = CONCAT1([:],d2_ch,bed)
	version_ch = version_ch.mix(d3_ch.version)

	col_ch = COLLECT_TO_FILE_01([suffix:".list"],d3_ch.vcf.collect())
	version_ch = version_ch.mix(col_ch.version)

	d4_ch = CONCAT2([:], col_ch.output.collect(),bed)
	version_ch = version_ch.mix(d4_ch.version)

	version_ch = MERGE_VERSION("concat vcfs", version_ch.collect())
emit:
	vcf = d4_ch.vcf
	index = d4_ch.index
	version = version_ch
}


