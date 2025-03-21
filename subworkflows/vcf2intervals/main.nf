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

include {VCF_TO_BED} from '../vcf2bed'
include {JVARKIT_VCF_TO_INTERVALS_01 as VCF_TO_INTERVALS} from '../../modules/jvarkit/vcf2intervals'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'



/**
 * call jvarkit vcf2intervals for one or more vcf
 *
 */
workflow JVARKIT_VCF_TO_INTERVALS_01 {
	take:
		vcf /* path to vcf, or file with .list suffix */
		bed /* limit to that BED or NO_FILE */
	main:	
		bed_ch = VCF_TO_BED(vcf)
	
		bed2ch = INTERSECT_BED(bed_ch.bed, bed)


		interval_vcf = bed2ch.bed.splitCsv(header:false,sep:'\t').
			map{T->[ T[0]+":"+((T[1] as int)+1)+"-"+T[2], file("NO_FILE") ,T[3],T[4]]}

		rch = VCF_TO_INTERVALS(interval_vcf)

		concat_ch = CONCAT_FILES_01(rch.bed.collect())

	emit:
		vcf2bed = bed_ch.bed /** original BED computed by VCF_TO_BED */
		bed = concat_ch.output /** sliced bed , this is what we want */
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


