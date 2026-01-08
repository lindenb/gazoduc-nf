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
include {isBlank          } from '../../../modules/utils/functions.nf'
include {verify           } from '../../../modules/utils/functions.nf'

process JVARKIT_VCF_CLUSTER {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id}"
	input:
		tuple val(meta  ),path(vcf)
	output:
		tuple val(meta ),path("*.vcf.gz"),optional:true,emit:vcf
		tuple val(meta ),path("*.vcf.gz.tbi"),optional:true,emit:tbi
		path("versions.yml"),emit:versions
	script:
	    def jvm = " -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
	    def jvarkit = "java ${jvm} -jar \${HOME}/jvarkit.jar" //TODO fix for conda
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		if(isBlank(args2)) throw new IllegalArgumentException("${task.process} args2 shouldn't be empty")
		def prefix =  task.ext.prefix?:"${meta.id}.bedcluster."
	"""
	mkdir -p TMP/OUT

	bcftools view  ${args1} "${vcf}" |\\
	${jvarkit} vcfcluster \\
		-o TMP/OUT \\
		--prefix "${prefix}" \\
	    ${args2} 
    
    mv TMP/OUT/*.vcf.gz ./ || true
    mv TMP/OUT/*.vcf.gz.tbi ./ || true

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
	"""


stub:
	def prefix = "${meta.prefix}.vcfcluster"
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz
"""
}
