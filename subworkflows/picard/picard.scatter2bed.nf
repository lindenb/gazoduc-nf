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
log.info("picard.scatter2bed.nf")
include { SCATTER_INTERVALS_BY_NS } from '../../modules/picard/picard.scatter.nf'
log.info("x2")
include { INTERVAL_LIST_TO_BED }  from '../../modules/picard/picard.interval2bed.nf' addParams(SORT:true, SCORE:500)
log.info("x3")

workflow SCATTER_TO_BED {
	take:
		genome
	main:
		version_ch = Channel.empty()
		scatter_ch = SCATTER_INTERVALS_BY_NS(genome)
		version_ch = version_ch.mix(scatter_ch.version)

		bed_ch = INTERVAL_LIST_TO_BED(scatter_ch.output)
	emit:
		bed = bed_ch.bed
		version = version_ch
	}
