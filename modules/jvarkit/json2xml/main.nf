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
process JSON_TO_XML {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
        tuple val(meta ),path(json)
	output:
		tuple val(meta), path("*.xml"),emit:xml
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
        def prefix = task.ext.prefix?:"${meta.id?:json.baseName}."
        def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
        def with_validation = (task.ext.with_validation?:true).toBoolean()
        def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released

 """
	hostname 1>&2
    mkdir -p TMP
	
    ${jvarkit} json2xml  \\
            ${args1} \\
            ${json}  > TMP/jeter.xml

    
    if ${with_validation}
    then
        xmllint  --nonet --noout  TMP/jeter.xml
    fi

    cp TMP/jeter.xml ${prefix}.xml

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
    xmllint: \$(xmllint --version 2>&1 |awk '(NR==1) {print \$NF;}')
EOF
	"""
	
stub: 
    def prefix = task.ext.prefix?:"${meta.id?:json.baseName}."
"""
touch versions.yml ${prefix}.xml
"""
}
