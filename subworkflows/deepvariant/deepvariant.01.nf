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
include {SAMTOOLS_SAMPLES} from '../samtools/samtools.samples.03.nf'
include {DEEP_VARIANT_CALL_01} from '../../modules/deepvariant/deepvariant.call.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {GLNEXUS_GENOTYPE_01} from '../../modules/glnexus/glnexus.genotype.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'

workflow DEEPVARIANT_01 {
	take:
		meta
		genomeId
		bams
		beds
	main:
		version_ch = Channel.empty()
		
		sn_ch = SAMTOOLS_SAMPLES([allow_multiple_references:false,with_header:true,allow_duplicate_samples:false],bams)
		version_ch = version_ch.mix(sn_ch.version)

		each_bed = beds.splitText().map{T->T.trim()}

		ch1 = sn_ch.rows.filter{T->T.genomeId.equals(genomeId)}
			.combine(each_bed)
			.map{T->T[0].plus("bed":file(T[1]))}
		
		ch2 = DEEP_VARIANT_CALL_01([:],ch1)
		version_ch = version_ch.mix(ch2.version)


		ch3 = ch2.output.map{T->T[0].plus("vcf":T[1],"gvcf":T[2])}.
			map{T->[T.bed,T.gvcf]}.
			groupTuple()

		ch4 = GLNEXUS_GENOTYPE_01([:],ch3)	
		version_ch = version_ch.mix(ch4.version)

		file_list_ch = COLLECT_TO_FILE_01([:],ch4.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(file_list_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output, file("NO_FILE"))
		version_ch = version_ch.mix(concat_ch.version)

                version_ch = MERGE_VERSION("Deep Variant",version_ch.collect())
			
	emit:
		version = version_ch
		vcf = concat_ch.vcf
	}

