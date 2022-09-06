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

include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {SAMTOOLS_STATS_01 as ST_STATS} from '../../modules/samtools/samtools.stats.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow SAMTOOLS_STATS_01 {
	take:
		meta
		reference
		references
		bams
	main:
		version_ch = Channel.empty()
		
		samples_ch = SAMTOOLS_SAMPLES_01(["with_header":"true"],reference,references,bams)
		version_ch = version_ch.mix(samples_ch.version)

		stats_ch = ST_STATS(meta,samples_ch.output.splitCsv(header:true,sep:"\t").map{T->T.plus("sample":T.new_sample)})
		version_ch = version_ch.mix(stats_ch.version)

                file_list_ch = COLLECT_TO_FILE_01([:], stats_ch.output.map{T->T[1]}.collect())
                version_ch = version_ch.mix(file_list_ch.version)

                multiqc_ch = MULTIQC_01(meta,file_list_ch.output.map{T->["files":T,"prefix":(meta.prefix?:"")]})
                version_ch = version_ch.mix(multiqc_ch.version)

		
		version_ch = MERGE_VERSION(meta, "ST stats", "samtools stats", version_ch.collect())
	emit:
		version = version_ch
		multiqc_zip = multiqc_ch.zip
	}

