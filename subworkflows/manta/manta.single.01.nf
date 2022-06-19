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
include { getKeyValue; getModules; getBoolean} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {MANTA_SINGLE_01} from  '../../modules/manta/manta.single.01.nf'
include {SURVIVOR_INSTALL}  from '../../modules/survivor/survivor.install.nf'

workflow MANTA_SINGLE_SV01 {
	take:
		meta
		reference
		bams
	main:
		version_ch = Channel.empty()

		all_samples_ch = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch= version_ch.mix(all_samples_ch.version)

		each_sample_bam_ch = all_samples_ch.out.splitCsv(header:false,sep:'\t')
	
		
		manta_ch = MANTA_SINGLE_01(meta.subMap(["prefix"]),reference,each_sample_bam_ch)
		version_ch= version_ch.mix(manta_ch.version.first())
	

		survivor_ch = SURVIVOR_INSTALL(meta)
		version_ch= version_ch.mix(survivor_ch.version)


		version_merge = MERGE_VERSION(meta,"manta","manta",version_ch.collect())
	emit:
		version = version_merge.version
	}

