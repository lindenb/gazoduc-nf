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
process TDB_LOADER {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(tdbloader_exe)
    tuple val(meta ),path(data_list)//one or more RDF data source
output:
    tuple val(meta ),path("TDB.*"),emit:datadir
	path("versions.yml"),emit:versions
script:
    def args1 = task.ext.args1?:""
    def jvm =  task.ext.jvm?:"-Djdk.xml.entityExpansionLimit=0 -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    def suffix = task.ext.suffix?:"${meta.id}"
    def data_array = (data_list instanceof List?data_list:[data_list])
"""
mkdir -p TMP
mkdir -p "TDB.${suffix}"

export JVM_ARGS="${jvm}"

${tdbloader_exe.toRealPath()} \\
   --loc=TDB.${suffix} \\
    ${data_array.join(" ")}

cat << EOF > versions.yml
"${task.process}"
	jena: TODO
EOF
"""

stub:
    def suffix = task.ext.suffix?:"${meta.id}"
"""
touch versions.yml
mkdir -p "TDB.${suffix}"
"""
}
