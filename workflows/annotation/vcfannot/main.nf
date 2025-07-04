/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../../subworkflows/jvarkit/jvarkit.vcf2intervals.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'


workflow {
	def refhash  = [ id: file(params.fasta) ]
	def fasta    = [ refhash, file(params.fasta) ]
	def fai      = [ refhash, file(params.fai) ]
	def dict     = [ refhash, file(params.dict) ]
	def pedigree = [ refhash, (params.pedigree ? file(params.pedigree) : []) ]
	def bed = [ refhash, (params.bed ? file(params.bed) : []) ]

	vcfs = Channel.fromPath(params.samplesheet)
		.splitCsv(header:true,sep=",")
		.map{
			if(it.id) return it;
			return it.plus(id:file(it.vcf).baseName);
			}
		.map { [id:it.id],file(params.vcf),file(params.index)]}

	ann_ch = ANNOTATE_VCF(
		[id:"annot"],
		fasta,
		fai,
		dict,
		vcfs
		bed,
		pedigree
		)
	}


workflow ANNOTATE_VCF {
	take:
		meta
		fasta
		fai
		dict
		vcf
		bed
		pedigree
	main:
		version_ch = Channel.empty()

		VCF_TO_BED(fasta,fai,dict,bed,vcf)

		ann_ch = ANNOTATE_VCF_01(meta,fasta,fai,dict, bed, ch1_ch)
		version_ch = version_ch.mix(ann_ch.versions)

		tofile_ch = COLLECT_TO_FILE_01([suffix:".list"], ann_ch.output.map{T->T.annot_vcf}.collect())
		version_ch = version_ch.mix(tofile_ch.versions)

		concat_ch = BCFTOOLS_CONCAT_01([:], tofile_ch.output, file("NO_FILE") )
		version_ch = version_ch.mix(concat_ch.version)

	emit:
		versions
		vcf = concat_ch.vcf
	}

runOnComplete(workflow);

process VCF_TO_BED {
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(bed)
	tuple val(meta ),path(vcf),path(idx)
output:
	tuple path("BEDS/*.bed"),val(meta),path(vcf),path(idx)
script:
	def method = task.ext.method?:""
	def has_bed= bed?true:false
"""
mkdir -p TMP
bcftools index -s "${vcf}" |\\
	awk -F '\t' '{printft("%s\t0\t%s\\n",\$1,\$2);}' |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter1.bed

if ${has_bed}
then
	${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\\
		grep -vE '^(browser|track)' |\\
		cut -f1,2,3 |\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed
	
	bedtools intersect -a TMP/jeter1.bed -b TMP/jeter2.bed |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter3.bed
	
	mv TMP/jeter3.bed TMP/jeter1.bed
fi

mkdir -p BEDS
jvarkit bedcluster ${method} -R ${fasta} -o BEDS TMP/jeter1.bed
"""
}

