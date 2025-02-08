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
include {BCFTOOLS_CONCAT_ALL} from '../bcftools/concat.all'


workflow TRUVARI_01 {
     take:
       	genome /* genome id */
        vcfs /* file containing the path to VCF files . One per line */
     main:

	vcf2bed_ch = VCF_TO_BED(vcfs)

	ch2 = vcfs.flatMap{[it,file(it.toRealPath()+(it.name.endsWith(".bcf")?".csi":".tbi"))]}.collect()	

	perctg_ch = PER_CONTIG(genome, vcf2bed_ch.chromosomes.splitText().map{it.trim()},ch2)

	concat_ch = BCFTOOLS_CONCAT_ALL(perctg_ch.output , file("NO_FILE"))

	emit:
		output = concat_ch.output
	}

process PER_CONTIG {
    tag "${contig}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../conda/truvari.01.yml"
    input:
	path(genome)
	val(contig)
	path("VCFS/*")
    output:
	tuple path("truvari.${contig}.bcf"),path("truvari.${contig}.bcf.csi"),emit: output
    script:
	def reference = genome.find{it.name.endsWith("a")}
	def extraCmd = task.ext.args1?:""
	def extraBcftools = task.ext.args2?:""
    """
	hostname 1>&2
	mkdir -p TMP

	find VCFS -type l \\( -name "*.vcf.gz" -o -name "*.bcf" \\) |\\
		while read F
		do
			bctools query -l "\${F}" | head -n 1 | tr "\\n" "\t" >> TMP/jeter.txt
			echo "\${F}" >> TMP/jeter.txt
		done
	
	sort -T TMP -t '\t' -k1,1 TMP/jeter.txt | cut -f 2 > TMP/jeter.list

	# bug FORMAT pour dragen
	bcftools merge --force-samples --filter-logic '+' --regions "${contig}" --file-list TMP/jeter.list  -m none -O u |\
		${extraBcftools.isEmpty()?"":"bcftools view -O u ${extraBcftools} |"} \
		bcftools annotate --force -x 'FORMAT/SR' -O z -o TMP/merged.vcf.gz
	bcftools index --tbi TMP/merged.vcf.gz

	# invoke truvari
	truvari collapse ${extraCmd} --reference "${reference}" -i "TMP/merged.vcf.gz" -c TMP/collapsed.vcf.gz |\
		bcftools +fill-tags -O u -- -t AN,AC,AF |\
		bcftools sort -T TMP -O b -o truvari.${contig}.bcf

	bcftools index -f truvari.${contig}.bcf
    """
   }
