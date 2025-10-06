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
process BED_TO_INTERVAL_LIST {
tag "${meta.id?:interval_list.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
when:
    task.ext.when == null || task.ext.when
input:
	tuple val(meta1),path(dict)
	tuple val(meta),path(bed)
output:
    tuple val(meta),path("*.interval_list"),emit:interval_list
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:bed.baseName
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir -p TMP
gatk --java-options "${jvm}" BedToIntervalList \\
    --INPUT "${bed}" \\
    --SD "${dict}" \\
    --O ${prefix}.interval_list


cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""

stub:
"""
touch versions.yml ${meta.id}.interval_list
"""
}
