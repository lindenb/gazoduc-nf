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
process BED_STATS {
	label "process_single"
	tag "${meta.id?:""} "
	afterScript "rm -rf TMP"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	input:
        tuple val(meta ),path("BED/*")
	output:
		tuple val(meta), path("*mqc.json",arity:"0..*"),emit:multiqc
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
        def jvm =  task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
        def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"// TODO update when conda released
		def prefix = task.ext.prefix?:"${meta.id}"
 """
    hostname 1>&2
    mkdir -p TMP
	find BED \\( -name "*.bed" -o -name "*.bed.gz" \\) > TMP/jeter.list

    ${jvarkit} bedstats  \\
            ${args1} \\
           -o TMP/STATS \\
           TMP/jeter.list 

mv -v TMP/STATS/*.json ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""

stub: 
    def prefix = task.ext.prefix?:"${meta.id}.dict2bed"
"""
touch versions.yml 
"""
}
