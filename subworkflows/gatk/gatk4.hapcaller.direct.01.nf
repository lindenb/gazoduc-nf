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
include {LINUX_SPLIT as LINUX_SPLIT1; LINUX_SPLIT as LINUX_SPLIT2} from '../../modules/utils/split.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {GATK4_HAPCALLER_DIRECT} from '../../modules/gatk/gatk4.hapcaller.direct.01.nf'
include {SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'
include {BCFTOOLS_MERGE_BED_01} from '../../modules/bcftools/bcftools.merge.bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {ANNOTATE_VCF_01} from '../annotation/annotation.vcf.01.nf'

workflow GATK4_HAPCALLER_DIRECT_01 {
	take:
		meta
		genomeId
		bams
		beds
	main:

		version_ch = Channel.empty()

		if(beds.name.equals("NO_FILE"))
			{
			log.warn("no beds defined. splitting the genome into parts using picard scatterreference")
			acgt_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"100"],params.genomes[genomeId].fasta)
			version_ch = version_ch.mix(acgt_ch.version)

			acgt_split_ch = LINUX_SPLIT1(["method":"--lines=1","suffix":".bed"],acgt_ch.bed)
			version_ch = version_ch.mix(acgt_split_ch.version)

			capture_bed_ch = acgt_split_ch.output.splitText().map{it.trim()}
			}
		else
			{
			capture_bed_ch = Channel.fromPath(beds).splitText().map{it.trim()}
			}

		version_ch = Channel.empty()
		split_bams_ch = LINUX_SPLIT2(["method":"--lines=${params.nbams?:"10"}","suffix":".list"],bams)
		version_ch = version_ch.mix(split_bams_ch.version)
		
		each_bam_list = split_bams_ch.output.splitText().map{it.trim()}

		row_ch = capture_bed_ch.combine(each_bam_list).map{T->[
					"bed":T[0],
					"bams":T[1],
					"genomeId":genomeId,
					"pedigree":params.pedigree?:"",
					"mapq":(params.mapq?:-1),
					"extraHC":(params.extraHC?:"")
					]}

		callbed_ch = GATK4_HAPCALLER_DIRECT([:],row_ch)
		version_ch = version_ch.mix(callbed_ch.version.first())

		bed2vcf_ch  = callbed_ch.bedvcf.map{T->[T[0].bed,T[1]]}.groupTuple()

		bedvcf_ch = BCFTOOLS_MERGE_BED_01([:], bed2vcf_ch)
		version_ch = version_ch.mix(bedvcf_ch.version)

		if(params.with_annotations) {
			ann_ch = ANNOTATE_VCF_01([:], genomeId, bedvcf_ch.bedvcf.map{T->[vcf:T[1],bed:T[0]]})
			version_ch = version_ch.mix(ann_ch.version)		
			vcfs_ch = ann_ch.output.map{T->T.annot_vcf}
			}
		else
			{
			vcfs_ch = bedvcf_ch.bedvcf.map{T->T[1]}
			}


		cut_ch = COLLECT_TO_FILE_01([:], vcfs_ch.collect())
		version_ch = version_ch.mix(cut_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:],cut_ch.output, file("NO_FILE"))
		version_ch = version_ch.mix(concat_ch.version)

		version_ch = MERGE_VERSION("gatk4Direct",version_ch.collect())

	emit:
		version = version_ch
		vcf = concat_ch.vcf
		index = concat_ch.index
	}

