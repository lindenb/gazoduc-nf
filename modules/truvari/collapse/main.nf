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
process TRUVARI_COLLAPSE {
    label "process_short"
	tag "${meta.id}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/truvari.01.yml"
    input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta ),path("VCFS/*")
    output:
		tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit: vcf
		path("versions.yml"),emit:versions
    script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		def args3 = task.ext.args3?:""
		def awk_expr = task.ext.awk_expr?:"^(chr)?[0-9XY]+\$"
		def prefix= task.ext.prefix?:meta.id
		def merge_mode = task.ext.merge_mod?:"id"
    """
	hostname 1>&2
	mkdir -p TMP
	set -x
	find VCFS -type l \\( -name "*.vcf.gz" -o -name "*.bcf" \\) |\\
		while read F
		do
			bcftools query -l "\${F}" | head -n 1 | tr "\\n" "\t" >> TMP/jeter.txt
			echo "\${F}" >> TMP/jeter.txt
		done
	
	sort -T TMP -t '\t' -k1,1 TMP/jeter.txt | cut -f 2 > TMP/jeter.list

	awk -F '\t' '(\$1 ~ /${awk_expr}/ ) {printf("%s\t0\t%s\\n",\$1,\$2);}' "${fai}" |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter.bed

	bcftools merge \\
		--regions-file TMP/jeter.bed \\
		--threads ${task.cpus} \\
		${args1} \\
		--force-samples \\
		--filter-logic '+' \\
		--file-list TMP/jeter.list \\
		-m ${merge_mode} \\
		-O u \\
		-o TMP/jeter2.bcf

	mv TMP/jeter2.bcf TMP/jeter.bcf

	#  Truvari: Cannot compare multi-allelic records. Please split
	bcftools norm \\
		--threads ${task.cpus} \\
		--fasta-ref ${fasta} \\
		--multiallelics -any \\
		-O u \\
		-o TMP/jeter2.bcf \\
		TMP/jeter.bcf

	mv TMP/jeter2.bcf TMP/jeter.bcf

	# optional filter ?
	bcftools view --threads ${task.cpus} ${args2}  -O u -o TMP/jeter2.bcf TMP/jeter.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf

	# bug DRAGEN
	bcftools annotate --threads ${task.cpus} --force -x 'FORMAT/SR' -O z -o TMP/merged.vcf.gz TMP/jeter.bcf
	bcftools index --threads ${task.cpus}  --tbi TMP/merged.vcf.gz

	# invoke truvari
	truvari collapse ${args3} --reference "${fasta}" -i "TMP/merged.vcf.gz" -c TMP/collapsed.vcf.gz |\\
		bcftools view -O u -o TMP/jeter.bcf
	
	bcftools +fill-tags --threads ${task.cpus}  -O u -o TMP/jeter2.bcf TMP/jeter.bcf -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
	mv TMP/jeter2.bcf TMP/jeter.bcf
		
	bcftools sort -T TMP/sort -O z -o TMP/${prefix}.vcf.gz TMP/jeter.bcf

	bcftools index -t --threads ${task.cpus} -f TMP/${prefix}.vcf.gz
	mv TMP/${prefix}.vcf.gz ./
	mv TMP/${prefix}.vcf.gz.tbi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	truvari: \$(truvari version)
END_VERSIONS
    """

stub:
def prefix= task.ext.prefix?:meta.id
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}
