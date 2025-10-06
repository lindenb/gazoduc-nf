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
process SCATTER_INTERVALS_BY_NS {
tag "${meta1.id?:fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
when:
    task.ext.when == null || task.ext.when

input:
	tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
    tuple val(meta1),path("*.interval_list"),emit:interval_list
    path("versions.yml"),emit:versions
script:
    def max_to_merge = task.ext.max_to_merge?:""
    def output_type = task.ext.output_type?:""
    if(max_to_merge.toString().trim().isEmpty()) throw new IllegalArgumentException("undefined max_to_merge for ${task.process}");
    if(output_type.toString().trim().isEmpty()) throw new IllegalArgumentException("undefined output_type for ${task.process}");
    def prefix = task.ext.prefix?:fasta.baseName
"""
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ScatterIntervalsByNs \\
    --REFERENCE "${fasta}" \\
    --MAX_TO_MERGE ${max_to_merge} \\
    --OUTPUT "TMP/jeter.interval_list" \\
    --OUTPUT_TYPE ${output_type}

mv "TMP/jeter.interval_list" ./${prefix}.interval_list

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""


stub:
"""
touch versions.yml "${meta1.id}.interval_list"
"""
}
