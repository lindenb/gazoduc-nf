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
include {moduleLoad;getModules} from '../../modules/utils/functions.nf'

/**
 * filter-out bed without variants in a VCF file
 */
process BED_WITH_VARIANTS_01 {
tag "${vcf.name} ${bed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	tuple path(vcf),path(bed)
output:
	tuple path(vcf),path(bed),path("bed_with_variants.bed"),path("excluded.bed"),emit:output
	path("version.xml"),emit:version
script:
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("bcftools bedtools")}
	set -x
	mkdir TMP	
	touch TMP/jeter.out.bed

	# prepare a BED to VCF

	## VCF is a list
	if ${vcf.name.endsWith(".list")} ; then
		cat "${vcf}" | while read -r V
		do
			bcftools query -f '%CHROM\t%POS0\t%END\\n' "${V}" | \
				sort -T TMP -t '\t' -k1,1 -k2,2n |\
				bedtools merge >> TMP/tmp1.bed
		done

		sort -T TMP -t '\t' -k1,1 -k2,2n TMP/tmp1.bed | bedtools merge > TMP/vcf.bed
	else
		# regular vcf
		bcftools query -f '%CHROM\t%POS0\t%END\\n'  '${vcf}' |\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\
			bedtools merge > TMP/vcf.bed
	fi

	# empty vcf ?
	if ! test -s TMP/vcf.bed ; then
		echo -e "xxxxxxxxxxxxxxx\t0\t1\tno_vcf.vcf.gz" > TMP/vcf.bed
	fi
	
	${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" | while read -r L
	do
		echo "\${L}" | cut -f1,2,3 > TMP/jeter.bed
		
		# loop over all VCF intersecting that region
		bedtools intersect -u -a TMP/vcf.bed -b TMP/jeter.bed | cut -f4 | while read -r V
		do
			# find variant in the current VCF
			bcftools query --regions-file TMP/jeter.bed -f 'OK\\n' "\${V}" | head -n1 > TMP/got.txt

			if test -s TMP/got.txt ; then
				echo "\${L}" >> TMP/jeter.out.bed
				break
			fi
		done
	done

	sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.out.bed > TMP/jeter.out2.bed

	${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" | sort -T TMP > TMP/a
	sort -T TMP TMP/jeter.out2.bed > TMP/b
	comm -13 TMP/a TMP/b > excluded.bed

	mv TMP/jeter.out2.bed bed_with_variants.bed
	


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs in multiple VCFs using bcftools</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="bed">${bed}</entry>
		<entry key="ctx.included">\$(wc -l < bed_with_variants.bed)</entry>
		<entry key="ctx.excluded">\$(wc -l < excluded.bed)</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
		<entry key="bedtools">\$( bedtools --version)</entry>
	</properties>
	EOF
	"""
	}
