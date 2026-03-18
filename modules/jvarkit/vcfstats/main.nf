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

process JVARKIT_VCFSTATS {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta1),path(sample2pop)
		tuple val(meta ),path(vcf),path(tbi)
	output:
		tuple val(meta ),path("*.json", arity: '0..*'),optional:true,emit:json
                tuple val(meta ),path("*.xml", arity: '0..*'),optional:true,emit:xml
		path("versions.yml"),emit:versions
	script:
        def args1 = task.ext.args1?:"-a"
        def args2 = task.ext.args2?:""
        def args3 = task.ext.args3?:""
        def prefix  = task.ext.prefix?:"${meta.id}.stats"
		def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
		def jvarkit = "java  ${jvm} -jar \${HOME}/jvarkit.jar" //todo fix me conda
	"""
	hostname 1>&2
    mkdir -p TMP/OUT

	${vcf instanceof List?"bcftools concat ${args1} -O u ${vcf.join(" ")} ":"cat ${vcf}"} |\\
        bcftools view ${args2} -O v |\\
		${jvarkit}  vcfstats \\
			${sample2pop?"--categories ${sample2pop}":""} \\
			${args3} \\
			-o TMP/OUT

    mv -v TMP/OUT/*.json ./ || true
    mv -v TMP/OUT/*.xml  ./ || true

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
EOF
	"""

stub:
 def prefix  = task.ext.prefix?:"${meta.id}.stats"
"""
touch versions.yml
"""

}
