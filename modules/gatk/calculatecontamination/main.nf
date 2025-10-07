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

process CALCULATE_CONTAMINATION {
tag "${meta.id?:""}"
label "process_single"
afterScript 'rm -rf TMP'
input:
    tuple val(meta),path(table1),path(opt_table2)
output:
    tuple val(meta),path("*.table"),emit:table
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

gatk --java-options "${jvm}" CalculateContamination \\
    -I "${table1}" \\
    ${opt_table2?"-matched \"${table2}\"":""} \\
    -O TMP/jeter.contamination.table

mv TMP/jeter.contamination.table "${meta.id}.contamination.table"

cat << EOF > version.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""
stub:
"""
touch versions.yml "${meta.id}.contamination.table"
"""
}
