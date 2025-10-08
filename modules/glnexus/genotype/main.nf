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



process GLNEXUS_GENOTYPE {
tag "${optional_bed?optional_bed.name:""} ${meta.id?:""}"
conda "${moduleDir}/../../../conda/glnexus.yml"
afterScript  "rm -rf TMP GLnexus.DB"
input:
	tuple val(meta1),path(optional_bed)
	tuple val(meta2),path(optional_config)
	tuple val(meta ),path("GVCF/*")
output:
	tuple val(meta ),path("*.bcf"),path("*.bcf.csi"),path(optional_bed),emit:vcf
	path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:""
	def config = task.ext.config?:""
	def prefix = task.ext.prefix?:(meta.id?:"glnexus")+(optional_bed?"."+optional_bed.baseName:"")+".glnexus"
	if(!optional_config && config.trim().isEmpty()) throw new IllegalArgumentException("${task.process} task.ext.config missing. eg.DeepVariantWGS")
"""
	hostname 1>&2
	mkdir TMP

	# glnexus doesn't like overlapping BED records
	if ${optional_bed?true:false}
	then
		cut -f1,2,3 "${optional_bed}" |\\
			sort -T TMP -t '\t' -k1,1 -k2,2n |\\
			bedtools merge > TMP/jeter.bed
	fi


	# check all samples in the same order
	find GVCF/ \\( -name "*.vcf.gz" -o -name ".gvcf.gz" \\) | while read F
	do
		bcftools query -l "\${F}" | head -n1 | tr "\\n" "," >> TMP/jeter.csv
		echo "\${F}" >> TMP/jeter.csv
	done
	
	# extract gvcf path
	sort -T TMP -k1,1 -t, TMP/jeter.csv |\\
		cut -d, -f2 > TMP/jeter.list
	test -s TMP/jeter.list
	
	# check no dup sample
	cut -f1 -d, TMP/jeter.csv | sort | uniq -d > TMP/jeter.dups
	cat TMP/jeter.dups 1>&2
	test ! -s TMP/jeter.dups
	
	
	glnexus_cli \\
	 	${args1} \\
		--config ${optional_config?"${optional_config}":"${config}"} \\
        ${optional_bed?"--bed TMP/jeter.bed":""} \\
		--threads ${task.cpus} \\
		--mem-gbytes ${task.memory.giga} \\
		--list TMP/jeter.list > TMP/jeter.vcf.gz
	
	bcftools +fill-tags ${args2} --threads ${task.cpus} -O b -o TMP/jeter.bcf TMP/jeter.vcf.gz -- -t  AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
	bcftools index --threads ${task.cpus} "TMP/jeter.bcf"


mv TMP/jeter.bcf ./${prefix}.bcf
mv TMP/jeter.bcf.csi ./${prefix}.bcf.csi

cat << EOF > versions.yml
"${task.process}":
    glnexus : "\$( bcftools view --header-only ${prefix}.bcf | grep -m1 -o GLnexusVersion=.*)"
EOF

"""
stub:
def prefix = task.ext.prefix?:(meta.id?:"glnexus")+(optional_bed?"."+optional_bed.baseName:"")+".glnexus"
"""
touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi
"""
}

