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
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(arq_exe)
    tuple val(meta2),path(data_list)//one or more RDF data source
    tuple val(meta ),path(query)
output:
    tuple val(meta ),path("*.arq.*"),emit:output
	path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def jvm =  task.ext.jvm?:"-Djdk.xml.entityExpansionLimit=0 -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def prefix = task.ext.prefix?:"${meta.id}"
    def suffix = task.ext.suffix?:"txt"
    def cmd1 = task.ext.cmd1?:"cat"
    def cmd2 = task.ext.cmd2?:"cat"
    def data_array = (data_list instanceof List?data_list:[data_list])
"""
mkdir -p TMP

#give a chance to alter the query, with seq or whatever
cat "${query}" | ${cmd1} > TMP/query.sparql

export JVM_ARGS="${jvm}"

${arq_exe.toRealPath()} \\
    ${args1} \\
    ${data_array.collect{"--data=${it}"}.join(" ")} \\
    --query=TMP/query.sparql | ${cmd2} > "TMP/jeter.tmp"

if ${suffix.endsWith(".gz")}
then
    gzip "TMP/jeter.tmp"
    mv "TMP/jeter.tmp.gz" "${prefix}.arq.${suffix}"
else
    mv "TMP/jeter.tmp" "${prefix}.arq.${suffix}"
fi 

cat << EOF > versions.yml
"${task.process}"
	jena: TODO
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
    def suffix = task.ext.suffix?:"txt"
"""
touch versions.yml "${prefix}.arq.${suffix}"
"""
}
