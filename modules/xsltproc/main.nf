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
process XSLTPROC {
	label "process_single"
	tag "${meta.id?:""} ${stylesheet.name}"
	afterScript "rm -rf TMP"
	input:
        tuple val(meta1),path(stylesheet)
        tuple val(meta ),path(xml)
	output:
		tuple val(meta), path("*.xslt.*"),emit:xml
		path("versions.yml"),emit:versions
	script:
	def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
        def prefix = task.ext.prefix?:"${meta.id?:xml.baseName}.${meta1.id?:stylesheet.baseName}"
        def suffix = task.ext.suffix?:"xml"
 """
	hostname 1>&2
    mkdir -p TMP
	
    xsltproc ${args1} ${args2} -o TMP/jeter.out "${stylesheet.toRealPath()}" "${xml}"

    if ${suffix.endsWith(".gz")}
    then
        gzip --best < TMP/jeter.out > "${prefix}.xslt.${suffix}"
    else
        mv TMP/jeter.out "${prefix}.xslt.${suffix}"
    fi

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
    xsltproc: \$(xsltproc --version 2>&1 |awk '(NR==1) {print \$NF;}')
EOF
	"""
	
stub: 
    def prefix = task.ext.prefix?:"${meta.id?:xml.baseName}.${meta1.id?:stylesheet.baseName}"
    def suffix = task.ext.suffix?:"xml"
"""
touch versions.yml ${prefix}.xslt.${suffix}
"""
}
