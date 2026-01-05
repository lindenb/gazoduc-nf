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

process JVARKIT_VCF_SET_DICTIONARY {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta1),path(dict_source)
		tuple val(meta ),path(vcf),path(idx)
	output:
		tuple val(meta ),path("*.{bcf,vcf.gz}"),path("*.{bcf.csi,vcf.gz.tbi}"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def size  = task.ext.size?:10
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:"-n SKIP"
		def args3 = task.ext.args3?:""
		def prefix  = task.ext.prefix?:vcf.baseName+".setdict"
		def sort = (task.ext.sort?:false).toBoolean()
		def suffix = task.ext.extension?:".bcf"
		def is_bcf = suffix.endsWith("bcf")
	"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

	bcftools view  ${args1} "${vcf}" |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfsetdict \\
			${args2} \\
			--reference '${dict_source}' |\\
	${sort?"bcftools sort -T TMP/sort -O u --max-mem \"${task.memory.giga}G\"  |":""} \\
	bcftools view ${args3} \\
		-O ${is_bcf?"b":"z"} \\
		-o TMP/${prefix}.${is_bcf?"bcf":"vcf.gz"}
	

	bcftools index \\
		-f \\
		${is_bcf?"":"-t"} \\
		--threads ${task.cpus} \\
		TMP/${prefix}.${is_bcf?"bcf":"vcf.gz"}
	

	mv TMP/${prefix}.${is_bcf?"bcf":"vcf.gz"} ./
	mv TMP/${prefix}.${is_bcf?"bcf.csi":"vcf.gz.tbi"} ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
	"""

stub:
 def prefix  = task.ext.prefix?:"${vcf.baseName}.setdict"
"""
touch versions.yml  ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""

}
