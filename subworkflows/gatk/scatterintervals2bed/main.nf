/*

Copyright (c) 2026 Pierre Lindenbaum

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

include { SCATTER_INTERVALS_BY_NS  } from '../../../modules/gatk/scatter.intervals/main.nf'
include { INTERVAL_LIST_TO_BED }  from '../../../modules/gatk/interval.list2bed/main.nf'

workflow SCATTER_TO_BED {
	take:
		meta
		fasta
		fai
		dict
	main:
		versions = Channel.empty()
		
		SCATTER_INTERVALS_BY_NS(fasta,fai,dict)
		versions = versions.mix(SCATTER_INTERVALS_BY_NS.out.versions)

		INTERVAL_LIST_TO_BED(SCATTER_INTERVALS_BY_NS.out.interval_list)
		versions = versions.mix(INTERVAL_LIST_TO_BED.out.versions)
	emit:
		interval_list = SCATTER_INTERVALS_BY_NS.out.interval_list
		bed = INTERVAL_LIST_TO_BED.out.bed
		versions
	}
