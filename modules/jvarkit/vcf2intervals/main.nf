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

process JVARKIT_VCF_TO_INTERVALS {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id}"
	input:
		tuple val(meta),path(vcf),path(idx),path(optional_bed)
	output:
		tuple val(meta),path("*.bed"),emit:bed /* chrom/start/end/vcf/idx */
		path("versions.yml"),emit:versions
	script:
		def cpu2 = java.lang.Math.max(1,task.cpus-3)
		def prefix = task.ext.prefix?:"${meta.id}"
		def awk_expr = task.ext.awk_expr?:"(1==1)"
		def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP -Xmx${task.memory.giga}G"
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		
	"""
	hostname 1>&2
	mkdir -p TMP

	bcftools view -G --threads ${cpu2} ${args1} ${optional_bed?"--regions-file \"${optional_bed.name}\"":""} "${vcf}" |\\
	java ${jvm}  -jar ${HOME}/jvarkit.jar  vcf2intervals  --bed  ${args2} |\\
		awk  -F '\t' '${awk_expr} {printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,"${vcf.toRealPath()}","${idx?idx.toRealPath().toString():""}");}' > TMP/intervals.bed

	mv TMP/intervals.bed "${prefix}.bed"

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
	"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.bed"
"""
}
