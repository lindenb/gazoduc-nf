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
include {isBlank} from '../../modules/utils/functions.nf'
include {BCFTOOL_STATS_01} from '../../modules/bcftools/bcftools.stats.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

workflow BCFTOOLS_VCFS_STATS_01 {
	take:
		meta
		reference
		vcfs
	main:
		version_ch = Channel.empty()
		each_vcf= vcfs.splitText().map{file(it.trim())}

		vcfstat_ch = BCFTOOL_STATS_01(meta,reference,each_vcf)
		version_ch = version_ch.mix(vcfstat_ch.version)

		version2_ch = MERGE_VERSION(meta, "VCF Stats", "VCF Stats using bcftools", version_ch.collect())
		html = VERSION_TO_HTML(params,version2_ch)

		zip_ch = SIMPLE_ZIP_01(["compression_level":"0"],vcfstat_ch.output.map{T->T[1]}.mix(version2_ch).mix(html.html).collect())
	emit:
		zip = zip_ch.zip
		version = version2_ch		
	}
