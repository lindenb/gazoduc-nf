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

include {VCF_TO_BED               } from '../../modules/bcftools/vcf2bed'
include {JVARKIT_VCF_TO_INTERVALS } from '../../modules/jvarkit/vcf2intervals'
include {CONCAT_FILES_01          } from '../../modules/utils/concat.files.nf'
include {BEDTOOLS_INTERSECT       } from '../../modules/bedtools/intersect'



/**
 * call jvarkit vcf2intervals for one or more vcf
 *
 */
workflow VCF_TO_INTERVALS {
	take:
		meta
		vcf /* meta,vcf,idx*/
		bed /* meta,bed */
	main:
		versions = Channel.empty()
		bed_ch = VCF_TO_BED(vcf)
		versions = versions.mix(VCF_TO_BED.out.versions)
	
		ch1 = VCF_TO_BED.out.bed.
			map{it instanceof List?it:[it]}.
			flatMap{v->{
				def L1=[];
				for(X in v[1]) {
					L1.add([v[0],X]);
					}
				return L1;
			}}.
			combine(bed).
			map{[it[0],it[1],it[3]]}.
			view()

		BEDTOOLS_INTERSECT(
			[ [id:"nofai"],[]],
			ch1
			)
		versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)
		

		ch2 =vcf.join(BEDTOOLS_INTERSECT.output.bed).view()

		JVARKIT_VCF_TO_INTERVALS() 

		/*
		bed2ch = INTERSECT_BED(bed_ch.bed, bed)


		interval_vcf = bed2ch.bed.splitCsv(header:false,sep:'\t').
			map{T->[ T[0]+":"+((T[1] as int)+1)+"-"+T[2], file("NO_FILE") ,T[3],T[4]]}

		rch = VCF_TO_INTERVALS(interval_vcf)

		concat_ch = CONCAT_FILES_01(rch.bed.collect())
		*/
	emit:
		versions
		//vcf2bed = bed_ch.bed /** original BED computed by VCF_TO_BED */
		//bed = concat_ch.output /** sliced bed , this is what we want */
	}


process INTERSECT_BED {
tag "${vcfbed.name} / ${userbed.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(vcfbed)
	path(userbed)
output:
	path("intersect.bed"),emit:bed
script:
"""
hostname 1>&2
set -o pipefail

mkdir TMP

if ${userbed.name.equals("NO_FILE")} ; then
	cp "${vcfbed}" intersect.bed
else

	sort -t '\t' -T TMP -k1,1 -k2,2n "${vcfbed}" >  TMP/a.bed
	sort -t '\t' -T TMP -k1,1 -k2,2n "${userbed}" | cut -f1,2,3 | bedtools merge >  TMP/b.bed
	bedtools intersect -a TMP/a.bed -b TMP/b.bed |\
		sort -t '\t' -T TMP -k1,1 -k2,2n > intersect.bed
fi
"""
}


