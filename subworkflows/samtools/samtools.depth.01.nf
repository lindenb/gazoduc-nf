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

include {SAMTOOLS_SAMPLES02} from '../../modules/samtools/samtools.samples.02.nf'
include {assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_DEPTH_01 as ST_DEPTH} from '../../modules/samtools/samtools.depth.01.nf'

workflow SAMTOOLS_DEPTH_01 {
	take:
		meta
		reference
		references
		bams
		bed
	main:
		version_ch = Channel.empty()
		
		samples_ch = SAMTOOLS_SAMPLES02(
			meta.plus(["with_header":"true"]),
			["reference":reference,"references":references,"bams":bams]
			)
		version_ch = version_ch.mix(samples_ch.version)

		dp_ch = ST_DEPTH(meta,samples_ch.output.splitCsv(header:true,sep:"\t").
			map{T->T.plus("bed":bed)}
			)
		version_ch = version_ch.mix(dp_ch.version)

                concat_ch = CONCAT_FILES_01(["downstream_cmd":" | LC_ALL=C sort -T . | uniq "],dp_ch.output.map{T->T[1]}.collect())
                version_ch = version_ch.mix(concat_ch.version)
		
		version_ch = MERGE_VERSION(meta, "Samtools depth", "samtools depth", version_ch.collect())
	emit:
		version = version_ch
		output = concat_ch.output
	}

