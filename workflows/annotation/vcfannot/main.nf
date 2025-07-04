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
nextflow.enable.dsl=2

include {moduleLoad;runOnComplete} from '../../../modules/utils/functions.nf'
include {ANNOTATE_VCF_01} from '../../../subworkflows/annotation/annotation.vcf.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../../subworkflows/jvarkit/vcf2intervals/main.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'


workflow {
	version_ch = Channel.empty()
	def ref_hash = [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName
		]
	def fasta = [ ref_hash, file(params.fasta) ]
	def fai =   [ ref_hash, file(params.fai) ]
	def dict =  [ ref_hash, file(params.dict) ]
	def dict =  [ [id:"pedigree"], params.pedigree==null? [] : file(params.pedigree) ]

	

	vcf2intervals_ch = JVARKIT_VCF_TO_INTERVALS_01([distance: params.vcf2interval_distance,min_distance:200], vcf,file("NO_FILE"))
			version_ch = version_ch.mix(vcf2intervals_ch.version)

			ch1_ch = vcf2intervals_ch.bed.splitCsv(sep:'\t',header:false).map{T->[
				vcf: file(T[3]),
				vcf_index: file(T[3]+ (T[3].endsWith(".vcf.gz")?".tbi":".csi")),
				bed: file("NO_FILE"),
				interval : T[0]+":"+((T[1] as int)+1)+"-"+T[2]
				]}

		ann_ch = ANNOTATE_VCF_01(genomeId, bed, ch1_ch)
		version_ch = version_ch.mix(ann_ch.version)

		tofile_ch = COLLECT_TO_FILE_01([suffix:".list"], ann_ch.output.map{T->T.annot_vcf}.collect())
		version_ch = version_ch.mix(tofile_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:], tofile_ch.output, file("NO_FILE") )
		version_ch = version_ch.mix(concat_ch.version)

		version_ch = MERGE_VERSION("vcfannot",version_ch.collect())
	emit:
		version = version_ch
		output = concat_ch.vcf
	}

runOnComplete(workflow);



