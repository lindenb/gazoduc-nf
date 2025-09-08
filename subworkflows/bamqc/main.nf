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
include {MOSDEPTH        } from '../../modules/mosdepth'
include {SAMTOOLS_STATS  } from '../../modules/samtools/stats'

workflow BAM_QC {
take:
	meta
	fasta
	fai
	dict
	bed /* meta, bed */
	bams /* meta, bam, bai */
main:
	versions_ch = Channel.empty()
	reports_ch  = Channel.empty()

	MOSDEPTH(
		fasta,
		fai,
		bams.combine(bed).map{[it[0],it[1],it[2],it[4]]} ,
		)
	versions_ch = versions_ch.mix(MOSDEPTH.out.versions)


        mosdepth_global  = MOSDEPTH.out.global_txt
        mosdepth_summary = MOSDEPTH.out.summary_txt
        mosdepth_regions = MOSDEPTH.out.regions_txt

	reports_ch = reports_ch
		.mix(mosdepth_global)
		.mix(mosdepth_summary)
		.mix(mosdepth_regions)

    SAMTOOLS_STATS(fasta,fai,bed,bams)
	versions_ch = versions_ch.mix(SAMTOOLS_STATS.out.versions)
	reports_ch = reports_ch.mix(SAMTOOLS_STATS.out.stats)

emit:
	versions = versions_ch
	multiqc = reports_ch
	mosdepth_global
	mosdepth_summary
	mosdepth_regions
}


