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

include {SAMTOOLS_SAMPLES02} from '../../modules/samtools/samtools.samples.02.nf'
include {assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {SAMTOOLS_FLAGSTATS_01} from '../../modules/samtools/samtools.flagstats.01.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow BAM_FLAGSTATS_01 {
	take:
		meta
		reference
		references
		bams
	main:
		version_ch = Channel.empty()
		
		samples_ch = SAMTOOLS_SAMPLES02(["with_header":"true"],["reference":reference,"bams":bams,"references":references])
		version_ch = version_ch.mix(samples_ch.version)

		flags_ch = SAMTOOLS_FLAGSTATS_01(meta,samples_ch.output.splitCsv(header:true,sep:"\t"))
		version_ch = version_ch.mix(flags_ch.version)

                multiqc_ch = MULTIQC_01(meta,flags_ch.output.map{T->T[1]}.collect())
                version_ch = version_ch.mix(multiqc_ch.version)

		
		
		version_ch = MERGE_VERSION(meta, "Flagstats", "samtools flagstats", version_ch.collect())
	emit:
		version = version_ch
		multiqc_zip = multiqc_ch.zip
	}

