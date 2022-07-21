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
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {MOSDEPTH_DOWNLOAD_01} from '../../modules/mosdepth/mosdepth.downoad.01.nf'
include {MOSDEPTH_RUN_01} from '../../modules/mosdepth/mosdepth.run.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow MOSDEPTH_BAMS_01 {
	take:
		meta
		reference
		bams
	main:
		version_ch  = Channel.empty()
		bams_ch = SAMTOOLS_SAMPLES01(meta,reference,bams)
		version_ch = version_ch.mix(bams_ch.version)	

		ch1 = bams_ch.output.splitCsv(header:false,sep:'\t').map{T->[
			"reference": reference,
			"sample": T[0],
			"bam": T[1],
			"mapq": (meta.mapq?:1)
			]}
		mosdepth_ch = MOSDEPTH_DOWNLOAD_01(meta)
		version_ch = version_ch.mix(mosdepth_ch.version)	
	
		ch2 = MOSDEPTH_RUN_01(meta, mosdepth_ch.executable,ch1)
		version_ch = version_ch.mix(ch2.version)	

		ch3 = MERGE_VERSION(meta, "Mosdepth", "Mosdepth", version_ch.collect())
	emit:
		version = ch3
		summary = ch2.summary
		globaldist = ch2.globaldist
		perbase = ch2.perbase
	}

