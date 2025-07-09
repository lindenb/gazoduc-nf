/*

Copyright (c) 2024 Pierre Lindenbaum

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
process INTERVAL_LIST_TO_BED {
tag "${meta.id?:interval_list.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(interval_list)
output:
    tuple val(meta),path("*.bed"),emit:bed
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:interval_list.baseName
    def expr = task.ext.awk_filter?:""
"""
mkdir -p TMP
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" IntervalListToBed \\
    --INPUT "${interval_list}" \\
    --OUTPUT "TMP/jeter.bed" \\
    --SORT true

if ${!expr.trim().isEmpty()}
then
    awk -F '\t' '(${expr})' "TMP/jeter.bed" > "TMP/jeter2.bed"
    mv "TMP/jeter2.bed" "TMP/jeter.bed"
fi

mv TMP/jeter.bed ./${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""
}
