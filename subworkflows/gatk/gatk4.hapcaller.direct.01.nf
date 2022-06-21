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
import {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
import {GATK4_HAPCALLER_DIRECT} from '../../modules/gatk/gatk4.hapcaller.direct.01.nf'

workflow GATK4_HAPCALLER_DIRECT_01 {
	take:
		meta
		reference
		references
		bams
		beds
		dbsnp
		pedigree
	main:
		version_ch = Channel.empty()
		split_bams_ch = LINUX_SPLIT(["method":"--lines=${meta.nbams?:"10"}","suffix":".list"],bams)
		version_ch = version_ch.mix(split_bams_ch.version)
		
		each_bam_list = split_bams_ch.output.splitCsv(header: false,sep:'\t',strip:true).map{T->T[0]}
		capture_bed_ch = beds.splitCsv(header: false,sep:',',strip:true)
		row_ch = capture_bed_ch.combine(capture_bed_ch).
				map{T->[
					"bed":T[0],
					"reference":T[1],
					"bam":T[2],
					"dbsnp":dbsnp,
					"pedigree":pedigree,
					"mapq":(meta.mapq?:-1),
					"extraHC":(meta.extraHC?:"")
					]}

		callbed_ch = GATK4_HAPCALLER_DIRECT(meta,row_ch)
		version_ch = version_ch.mix(callbed_ch.version.first())

		ped2vcf_ch  = callbed_ch.bedvcf.groupTuple()
		mergedbed_ch = MERGE_BED(meta, reference, ped2vcf_ch)
		version_ch = version_ch.mix(mergedbed_ch.version)

		COLLECT_TO_FILE(meta,mergedbed_ch.collect())

		BCFTOOLS_CONCAT_01(meta,,"")

	emit:
		version = version_ch
	}

