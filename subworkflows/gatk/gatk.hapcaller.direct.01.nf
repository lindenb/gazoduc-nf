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
import {LINUX_SPLIT} from '../../modules/utils/split.nf'

workflow GATK_HAPCALLER_DIRECT_01 {
	take:
		meta
		reference
		references
		bams
		beds
	main:
		version_ch = Channel.empty()
		split_bams_ch = LINUX_SPLIT(["method":"--lines=${params.nbams?:"10"}","suffix":".list"],bams)
		version_ch = version_ch.mix(split_bams_ch.version)
		
		each_bam_list = split_bams_ch.output.splitCsv(header: false,sep:'\t',strip:true).map{T->T[0]}
		
	emit:
		output = merge_ch.out
		version = version_ch
	}

