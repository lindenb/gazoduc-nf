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

process JENA_ARQ {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(arq_exe)
    tuple val(meta2),path(data)
    tuple val(meta ),path(query)
output:
    tuple val(meta ),path("*.result.*"),emit:output
	path("versions.yml"),emit:versions
script:
    def version = task.ext.version?:"5.6.0"
    def jvm =  task.ext.jvm?:"-Djdk.xml.entityExpansionLimit=0 -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def prefix = task.ext.prefix?:"${meta.id}"
    def suffix = task.ext.suffix?:"txt"
    def cmd1 = task.ext.cmd1?:"cat"
"""
mkdir -p TMP
${arq_exe} ${jvm} --data=${data} --query=${query} | ${cmd1} > "${prefix}.result.${suffix}"

cat << EOF > versions.yml
"${task.process}"
	jena: TODO
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
    def suffix = task.ext.suffix?:"txt"
"""
touch versions.yml "${prefix}.result.${suffix}"
"""
}
