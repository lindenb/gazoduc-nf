/*

Copyright (c) 2025 Pierre Lindenbaum

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

process GATK4_GATHER_BQSR {
tag "${meta.id?:""}"
label "process_single"
afterScript 'rm -rf TMP'
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
	tuple val(meta),path("TABLES/*")
output:
    tuple val(meta),path("*.recal.table"),emit:table
    path("versions.yml"),emit:versions
script:
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

find TABLES/ -name "*.table" |\\
	awk '{printf("-I %s\\n",\$0);}' > TMP/arguments.list

gatk --java-options "${jvm}" GatherBQSRReports \\
	--arguments_file  TMP/arguments.list \\
	-O "${meta.id}.recal.table" 

cat << EOF > version.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""
stub:
"""
touch versions.yml "${meta.id}.recal.table" 
"""
}
