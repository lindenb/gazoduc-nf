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

process JVARKIT_VCFFILTERJDK {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id}"
	input:
		tuple val(meta1 ),path(script_file) //optional pr define expression
		tuple val(meta2 ),path(opt_pedigree)
		tuple val(meta  ),path(vcf)
	output:
		tuple val(meta ),path("*.vcf.gz"),optional:true,emit:vcf
		path("versions.yml"),emit:versions
	script:
        def has_script_file = (script_file?true:false)
		def expression = task.ext.expression?:""
        verify(
		(isBlank(expression) && has_script_file) ||
		(!has_script_file && !isBlank(expression)),
		"${task.process}: script file=${has_script_file} XOR expression=${expression} must be defined")
	def jvm = " -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
        def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
        def args3 = task.ext.args3?:""
        def prefix =  task.ext.prefix?:"${meta.id}.filterjdk"
	"""
	mkdir -p TMP

	bcftools view  ${args1} "${vcf}" |\\
	jvarkit ${jvm} vcffilterjdk \\
        ${args2} \\
        ${opt_pedigree?"--pedigree \"${opt_pedigree}\"":""} \\
        ${has_script_file?"--script \"${script_file}\"":""} \\
        ${isBlank(expression)?"":"--expression '${expression}'"}  | \\
	bcftools view ${args3} -O z -o TMP/jeter.vcf.gz

    if test \$(bcftools query -f '\\n'  TMP/jeter.vcf.gz |wc -l) -gt 0
	then
		mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
	fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
	"""


stub:
	def prefix = "${meta.prefix}"
"""
touch versions.yml ${prefix}.vcf.gz 
"""
}
