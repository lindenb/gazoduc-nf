/*

Copyright (c) 2025 Pierre Lindenbaum

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

include {BASE_RECALIBRATOR       } from '../../../modules/gatk/baserecalibrator'
include {APPLY_BQSR              } from '../../../modules/gatk/applybqsr'

workflow BQSR {
take:
	meta
	fasta
	fai
	dict
	known_variants // meta,vcf_files
	bams
main:
	versions = Channel.empty()

	BASE_RECALIBRATOR(
		fasta,
		fai,
		dict,
		known_variants,
		bams
		)
	versions = versions.mix( BASE_RECALIBRATOR.out.versions)
	
	APPLY_BQSR(
		fasta,
		fai,
		dict,
		BASE_RECALIBRATOR.out.table
		);
	versions = versions.mix( APPLY_BQSR.out.versions)
emit:
	versions
	bams = APPLY_BQSR.out.bam