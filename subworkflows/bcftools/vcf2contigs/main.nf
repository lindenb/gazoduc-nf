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
include { VCF_TO_BED            } from '../../../modules/bcftools/vcf2bed'

workflow VCF_TO_CONTIGS {
take:
    metadata
    vcfs //[meta,vcf,tbi]
main:
    versions = Channel.empty()
    
	VCF_TO_BED(vcfs)
	versions = versions.mix(VCF_TO_BED.out.versions)


    vcf = VCF_TO_BED.out.bed
		.splitCsv(header:false,sep:'\t')
		.map{meta,row->[meta,row[0]]}
		.combine(vcfs)
		.filter{meta1,_contig,meta2,_vcf,_tbi->meta1.id==meta2.id}
		.map{meta1,ctg,_meta2,vcffile,tbi->[meta1.plus(contig:ctg),vcffile,tbi]}
emit:
    vcf
    versions
}