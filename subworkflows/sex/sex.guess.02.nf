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
include {getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SEX_GUESS_01} from './sex.guess.01.nf'

/**
guess sex from list of bams
*/
workflow SEX_GUESS_02 {
	take:
		meta
		reference
		references
		bams
	main:
		version_ch = Channel.empty()

		bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":true,"allow_multiple_references":true,"allow_duplicate_samples":true]),reference,references,bams)
		version_ch = version_ch.mix(bams_ch.version)
		
		sex_ch= SEX_GUESS_01(meta,bams_ch.output)
		version_ch = version_ch.mix(sex_ch.version)

		version_ch = MERGE_VERSION(meta, "SexGuess", "Guess Sex from Bams", version_ch.collect())
	emit:
		/* version */
		version= version_ch
		pdf = sex_ch.pdf
		output = sex_ch.output
	}


