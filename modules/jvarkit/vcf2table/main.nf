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

process JVARKIT_VCF_TO_TABLE {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id}"
	input:
		tuple val(meta1),path(opt_pedigree)
		tuple val(meta ),path(vcf)
	output:
		tuple val(meta),path("*.{html,txt}"),emit:table
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
		def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP -Xmx${task.memory.giga}G"
        def prefix = task.ext.prefix?:"${meta.id}.table"
        def suffix = args2.contains("html") || "${task.process}".toLowerCase().contains("html")?"html":"txt"
"""
mkdir -p TMP

bcftools view ${args1} "${vcf}" |\\
jvarkit ${jvm} vcf2table  \\
    ${opt_pedigree?"--pedigree ${opt_pedigree}":""} \\
    ${args2} > TMP/jeter.table
    
mv -v TMP/jeter.table ${prefix}.${suffix}

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
EOF
	"""

stub:
    def args2 = task.ext.args2?:""
    def prefix = task.ext.prefix?:"${meta.id}"
    def suffix = args2.contains("html")?"html":"txt"
"""
touch versions.yml ${prefix}.${suffix}
"""
}
