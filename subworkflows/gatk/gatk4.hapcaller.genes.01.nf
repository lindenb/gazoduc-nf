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
include {GATK4_HAPCALLER_BED_01} from './gatk4.hapcaller.bed.01.nf'
include {BED_FOR_GENES_01} from '../../subworkflows/gtf/bed.for.genes.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

workflow GATK4_HAPCALLER_GENES_01 {
	take:
		meta
		reference
		references
		bams
		genelist
		dbsnp
		pedigree
	main:
		version_ch = Channel.empty()

		bed_genes_ch = BED_FOR_GENES_01(meta,reference,genelist)
		version_ch = version_ch.mix(bed_genes_ch.version)
		
		bed_genes_ch.bed.view{"BED for Genes is: $it ."}

		gatk_ch = GATK4_HAPCALLER_BED_01(meta,reference,references,bams,bed_genes_ch.bed,dbsnp,pedigree)
		version_ch = version_ch.mix(gatk_ch.version)

		version_ch = MERGE_VERSION(meta, "gatk4 gene", "call bams in a set of genes", version_ch.collect())
	emit:
                vcf = gatk_ch.vcf
                index = gatk_ch.index
		version = version_ch.version
	}

