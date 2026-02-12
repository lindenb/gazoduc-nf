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
include { VCF_TO_BED                   }  from '../../../modules/bcftools/vcf2bed'
include { BCFTOOLS_CONCAT_ALL          }  from '../../../subworkflows/bcftools/concat.all'

/** concat VCF per contig : IMPORTANT --regions meta.config must be set in processes under BCFTOOLS_CONCAT_ALL */
workflow BCFTOOLS_CONCAT_CONTIGS {
take:
	metadata
	vcfs // tuple [meta,vcf,idx]
main:
	versions = Channel.empty()
	multiqc = Channel.empty()

	vcfs = vcfs.map{meta,vcf,tbi->[meta.plus(id:vcf.toRealPath().toString().md5()),vcf,tbi]}

	VCF_TO_BED(vcfs)
	versions  = versions.mix(VCF_TO_BED.out.versions)

	ch1 = VCF_TO_BED.out.bed
		.splitCsv(header:false,sep:'\t')
		.map{meta,row->meta.plus(contig:row[0])}
		.combine(vcfs)
		.filter{meta1,meta2,_vcf,_tbi->meta1.id==meta2.id}
		.map{meta1,_meta2,vcf,tbi->[[id:meta1.contig+"${metadata.id?:""}",contig:meta1.contig],vcf,tbi]}

	BCFTOOLS_CONCAT_ALL(
		metadata,
		ch1
		)
	versions  = versions.mix(BCFTOOLS_CONCAT_ALL.out.versions)
	multiqc   = versions.mix(BCFTOOLS_CONCAT_ALL.out.multiqc)
	
emit:
	vcf = BCFTOOLS_CONCAT_ALL.out.vcf.view()
	versions
	multiqc
}

